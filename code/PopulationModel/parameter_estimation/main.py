#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
pcmp_unified.py

Unified Physiology-based Climate-driven Mosquito Population (PCMP) estimation
program for two species and a hierarchy of model configurations.

This single program reproduces the full model-variant series used in the
manuscript, from the baseline PCMP up to the full local-process model, for
either mosquito species, under either time-integration scheme, all estimated
with a common Differential Evolution (DE) optimizer.

------------------------------------------------------------------------------
What you choose at run time
------------------------------------------------------------------------------
1. Species:
       Ae  -> Aedes albopictus  (egg diapause)
       Cx  -> Culex pipiens     (adult diapause)
   The two species share an identical stage-structured framework, climate
   forcing, development/mortality functions, washout formulation, and
   estimation procedure. They differ ONLY in the diapause transition
   (which life stage overwinters) and in two parameter bounds.

2. Integration scheme:
       E -> explicit Euler
       R -> 4th-order Runge-Kutta (RK4)
   The governing equations are identical; only the numerical step differs.

3. Model configuration (which local-process components are active):

   Baseline
       B_PCMP            threshold (old) washout, s_obs fixed at 0.01,
                         mortality multipliers fixed at 1.0
   Single-component additions to the baseline
       M*W               + proportional (revised) washout
       M*S               + observation-scaling factor s_obs (estimated)
       M*M               + stage-specific mortality multipliers (estimated)
   Full model
       M*A               proportional washout + s_obs + mortality (all on)
   Leave-one-out from the full model (for contribution analysis)
       M*A-washout       full model with washout reverted to threshold type
       M*A-s_obs         full model with s_obs fixed
       M*A-mortality     full model with mortality fixed

   The "*" is the integration prefix: E (Euler) or R (RK4). For example,
   "MRA" is RK4 + all components (this corresponds to the manuscript's
   final working model), and "MEW" is Euler + washout only.

   Note on estimation: the original PCMP was historically estimated with an
   MCMC procedure. In this program EVERY configuration, including the
   baseline, is re-estimated under the SAME Differential Evolution optimizer
   so that all variants are compared on an identical maximum-likelihood
   footing. This is a deliberate methodological choice to keep the comparison
   fair; it is not a reproduction of the original MCMC estimation.

------------------------------------------------------------------------------
Output (per run)
------------------------------------------------------------------------------
A JSON file and a CSV file are written, in the same style as the original
single-configuration scripts, so that existing post-processing and figure
code can read them directly:

    <SPECIES>_<INTEGRATOR>_<CONFIG>_<POINT>_params.json
    output_paraest/<POINT>/<SPECIES_FULL>/<CONFIG_CODE>/..._est_week.csv

The JSON records the human-readable label, the full descriptive name, the
active-component flags, the raw and converted parameters, and the full set of
fit statistics, so the configuration is fully self-documenting.

------------------------------------------------------------------------------
Usage
------------------------------------------------------------------------------
Interactive:
    python pcmp_unified.py
    (you will be prompted for species, integrator, and configuration)

Command line:
    python pcmp_unified.py --species Ae --integrator R --config MA
        (--integrator R + --config MA  ->  label "MRA", the full RK4 model)
    python pcmp_unified.py --species Cx --integrator E --config B_PCMP
    python pcmp_unified.py --species Ae --integrator R --config all
        ("all" is NOT a model: it is a batch run mode that estimates every
         configuration (B_PCMP, MW, MS, MM, MA, and the three leave-one-out
         variants) in turn for the chosen species/integrator and writes a
         combined summary CSV. The full model MA is included exactly once,
         so "all" and "MA" are not duplicates: MA is one model, "all" is the
         batch that contains it.)

Quick test (small DE budget):
    python pcmp_unified.py --species Ae --integrator R --config MA \
        --de-maxiter 20 --de-popsize 8
"""

import os
import sys
import math
import json
import time
import argparse
from datetime import datetime

import numpy as np
import pandas as pd

try:
    from numba import njit
    NUMBA_AVAILABLE = True
except Exception:
    NUMBA_AVAILABLE = False

    def njit(*args, **kwargs):
        def deco(f):
            return f
        if args and callable(args[0]):
            return args[0]
        return deco

try:
    from scipy.optimize import differential_evolution
    from scipy.stats import pearsonr, spearmanr
    SCIPY_AVAILABLE = True
except Exception:
    SCIPY_AVAILABLE = False


# ============================================================
# GLOBAL SETTINGS (shared by both species)
# ============================================================

POINT = "tokyo"
OBSPOINT = "tokyo"

START_YEAR = 2003
YEAR_LENGTH = 11
WEEK_LENGTH = 52 * YEAR_LENGTH
WARMUP_START_YEAR = 1991
WARMUP_END_YEAR = 2002

# Sub-daily time step (shared by Euler and RK4 so that the two schemes are
# compared at the same temporal discretization).
DAY_LENGTH = 25
DELTA_T = 1.0 / DAY_LENGTH

# Differential Evolution settings (production defaults).
DE_MAXITER = 300
DE_POPSIZE = 15
DE_TOL = 1.0e-6
DE_SEED = 123
DE_POLISH = False
DE_WORKERS = 1
DE_UPDATING = "immediate"
DE_MUTATION = (0.5, 1.0)
DE_RECOMBINATION = 0.7

PRINT_EVERY_OBJECTIVE = 200
MAX_PENALTY = 1.0e12

# Fixed values used when a component is switched off.
FIXED_DIS_L = -2.0     # s_obs = 10**(-2) = 0.01  (original PCMP fixed value)
FIXED_M_MULT = 1.0     # mortality multipliers neutral
FIXED_PF_SCALAR_L = 0.0  # unused under threshold washout


# ============================================================
# PHYSIOLOGICAL CONSTANTS (shared)
# ============================================================

D0_TW_MIN = 3.14
D1_TW_MIN = 4.06
D2_TW_MIN = 5.31
D3_TA_MIN = 10.0

M1_TW_MAX = 35.04
M3_TA_MAX = 39.2
M0_TW_MIN = 14.64
M1_TW_MIN = 5.16
M2_TW_MIN = 5.16
M3_TA_MIN = 3.18

# Cx adult-diapause baseline mortality. For Aedes (egg diapause) the diapause
# mortality is tied to the egg mortality m0; for Culex (adult diapause) the
# overwintering adults instead carry a low fixed baseline mortality.
M4_BASE_CX = 0.005


# ============================================================
# PARAMETER DEFINITIONS
# ============================================================

ALL_PARAM_NAMES = [
    "p_b1",
    "p_b2",
    "p_ccl",
    "p_a1l",
    "p_a2l",
    "p_pf",
    "p_pf_scalar_l",
    "p_dis_l",
    "p_m_mult_egg",
    "p_m_mult_larva",
    "p_m_mult_pupa",
    "p_m_mult_adult",
    "p_m_mult_diapause",
]

# Bounds differ between species only for p_ccl and p_dis_l upper limits.
BOUNDS_AE = np.array([
    [10.0, 14.0], [10.0, 14.0], [0.0, 6.0], [-2.0, 10.0], [-2.0, 10.0],
    [0.5, 200.0], [-3.0, 10.0], [-5.0, 5.0],
    [0.01, 10.0], [0.01, 10.0], [0.01, 10.0], [0.01, 10.0], [0.01, 10.0],
], dtype=np.float64)

BOUNDS_CX = np.array([
    [10.0, 14.0], [10.0, 14.0], [0.0, 10.0], [-2.0, 10.0], [-2.0, 10.0],
    [0.5, 200.0], [-3.0, 10.0], [-5.0, 10.0],
    [0.01, 10.0], [0.01, 10.0], [0.01, 10.0], [0.01, 10.0], [0.01, 10.0],
], dtype=np.float64)

# Index groups (positions in ALL_PARAM_NAMES).
CORE_IDX = [0, 1, 2, 3, 4, 5]   # b1, b2, ccl, a1l, a2l, pf
PF_SCALAR_IDX = [6]             # proportional washout intensity
DIS_IDX = [7]                   # observation-scaling factor s_obs
MORT_IDX = [8, 9, 10, 11, 12]   # stage-specific mortality multipliers


# ============================================================
# SPECIES METADATA
# ============================================================

SPECIES_INFO = {
    "Ae": {
        "code": "Ae",
        "full": "Aedes_albopictus",
        "pretty": "Aedes albopictus",
        "diapause": "egg",
        "bounds": BOUNDS_AE,
        "obs_candidates": [
            "Aealbopictus_female", "Aealbo_female",
            "Aedes_albopictus_female", "Aedes_albopictus", "Aealbo",
        ],
    },
    "Cx": {
        "code": "Cx",
        "full": "Culex_pipiens",
        "pretty": "Culex pipiens",
        "diapause": "adult",
        "bounds": BOUNDS_CX,
        "obs_candidates": [
            "Cxpipiens_female", "Culex_pipiens_female",
            "Cx_pipiens_female", "Cxpipiens", "Culex_pipiens",
        ],
    },
}


# ============================================================
# MODEL CONFIGURATIONS
# ============================================================
# Each configuration is fully described by three boolean flags:
#   use_new_washout : proportional washout (True) vs threshold washout (False)
#   use_dis         : estimate s_obs (True) vs fixed at 0.01 (False)
#   use_mortality   : estimate mortality multipliers (True) vs fixed at 1.0
#
# The integration scheme (Euler/RK4) is orthogonal and is chosen separately.
# Labels below use a neutral "M" + integrator-letter convention; the actual
# label is built at run time with the integrator prefix substituted in.

CONFIG_FLAGS = {
    # key            (new_washout, dis,   mortality)
    "B_PCMP":        (False,       False, False),
    "MW":            (True,        False, False),  # washout only
    "MS":            (False,       True,  False),  # s_obs only
    "MM":            (False,       False, True),   # mortality only
    "MA":            (True,        True,  True),   # all on (full)
    "MA-washout":    (False,       True,  True),   # full minus washout
    "MA-s_obs":      (True,        False, True),   # full minus s_obs
    "MA-mortality":  (True,        True,  False),  # full minus mortality
}

CONFIG_ORDER = [
    "B_PCMP", "MW", "MS", "MM", "MA",
    "MA-washout", "MA-s_obs", "MA-mortality",
]

COMPONENT_LABELS = {
    "washout_new": "proportional washout",
    "washout_old": "threshold washout",
    "dis_on": "s_obs estimated",
    "dis_off": "s_obs fixed (0.01)",
    "mort_on": "mortality estimated",
    "mort_off": "mortality fixed (1.0)",
}


def build_label(config_key, integrator):
    """Human-readable short label, e.g. B_PCMP -> 'R_B_PCMP'? No: keep B_PCMP
    integrator-agnostic in name but record integrator separately. For the
    additive/full/LOO configs, prefix with M+integrator letter (E/R)."""
    if config_key == "B_PCMP":
        return f"B_PCMP[{integrator}]"
    # MW/MS/MM/MA/MA-... -> insert integrator letter after 'M'
    suffix = config_key[1:]  # drop leading 'M'
    return f"M{integrator}{suffix}"


def build_full_name(config_key, integrator, diapause):
    new_wash, dis_on, mort_on = CONFIG_FLAGS[config_key]
    integ = "Euler" if integrator == "E" else "RK4"
    parts = [f"{integ} integration"]
    parts.append(COMPONENT_LABELS["washout_new"] if new_wash else COMPONENT_LABELS["washout_old"])
    parts.append(COMPONENT_LABELS["dis_on"] if dis_on else COMPONENT_LABELS["dis_off"])
    parts.append(COMPONENT_LABELS["mort_on"] if mort_on else COMPONENT_LABELS["mort_off"])
    base = "PCMP baseline" if config_key == "B_PCMP" else "PCMP"
    return f"{base} ({diapause} diapause); " + "; ".join(parts) + "; Differential Evolution"


# ============================================================
# DATA LOADER
# ============================================================

def _read_climate_6(fname):
    df = pd.read_csv(fname, header=None)
    if df.shape[1] < 6:
        raise ValueError(f"Climate file has fewer than 6 columns: {fname}")
    df = df.iloc[:, :6]
    df.columns = ["TA", "TW", "DWK", "SMOI", "PREC", "ROFF"]
    return df


def load_input_arrays(species_full, obs_candidates, point=POINT, obs_point=OBSPOINT):
    ta_wup, tw_wup, dwk_wup, smoi_wup, prec_wup, roff_wup = [], [], [], [], [], []

    for yr in range(WARMUP_START_YEAR, WARMUP_END_YEAR + 1):
        fname = f"./../../input_data/{point}/{yr}/{point}{yr}_climdata2.csv"
        if os.path.exists(fname):
            df = _read_climate_6(fname)
            ta_wup.extend(df["TA"].values.tolist())
            tw_wup.extend(df["TW"].values.tolist())
            dwk_wup.extend(df["DWK"].values.tolist())
            smoi_wup.extend(df["SMOI"].values.tolist())
            prec_wup.extend(df["PREC"].values.tolist())
            roff_wup.extend(df["ROFF"].values.tolist())
        else:
            print(f"Warning: warmup climate file missing: {fname}")

    ta, tw, dwk, smoi, prec, roff = [], [], [], [], [], []
    data_year = [0] * YEAR_LENGTH

    for yr in range(START_YEAR, START_YEAR + YEAR_LENGTH):
        fname1 = f"./../../input_data/{point}/{yr}/{point}{yr}_climdata.csv"
        fname2 = f"./../../input_data/{point}/{yr}/{point}{yr}_climdata2.csv"
        if os.path.exists(fname1):
            fname = fname1
        elif os.path.exists(fname2):
            fname = fname2
        else:
            print(f"Warning: main climate file missing: {fname1} and {fname2}")
            continue
        df = _read_climate_6(fname)
        idx = yr - START_YEAR
        data_year[idx] = len(df)
        ta.extend(df["TA"].values.tolist())
        tw.extend(df["TW"].values.tolist())
        dwk.extend(df["DWK"].values.tolist())
        smoi.extend(df["SMOI"].values.tolist())
        prec.extend(df["PREC"].values.tolist())
        roff.extend(df["ROFF"].values.tolist())

    obs_adults_week = np.full(WEEK_LENGTH, -999.0, dtype=np.float64)
    species1 = None
    for s1 in obs_candidates:
        test_path = f"./../../measurement_data/{obs_point}/{species_full}/{obs_point}_{s1}_{START_YEAR}.tsv"
        if os.path.exists(test_path):
            species1 = s1
            print(f"Data file found: {test_path}")
            break
    if species1 is None:
        raise FileNotFoundError(
            f"No observation data found in ./../../measurement_data/{obs_point}/{species_full}/"
        )

    week_plus = 0
    for yr in range(START_YEAR, START_YEAR + YEAR_LENGTH):
        fname = f"./../../measurement_data/{obs_point}/{species_full}/{obs_point}_{species1}_{yr}.tsv"
        if os.path.exists(fname):
            df = pd.read_csv(fname, sep="\t", header=None, names=["week", "pop"], engine="python")
            if yr != START_YEAR:
                week_plus += 52
            for _, row in df.iterrows():
                wk = int(row["week"])
                pop = float(row["pop"])
                wc = wk + week_plus - 1
                if 0 <= wc < WEEK_LENGTH:
                    obs_adults_week[wc] = pop
        else:
            print(f"Warning: observation file missing: {fname}")
            if yr != START_YEAR:
                week_plus += 52

    to_np = lambda v: np.asarray(v, dtype=np.float64)
    return (
        to_np(ta_wup), to_np(tw_wup), to_np(dwk_wup), to_np(smoi_wup), to_np(prec_wup), to_np(roff_wup),
        to_np(ta), to_np(tw), to_np(dwk), to_np(smoi), to_np(prec), to_np(roff),
        obs_adults_week,
        np.asarray(data_year, dtype=np.int64),
    )


# ============================================================
# SHARED PHYSIOLOGICAL FUNCTIONS (njit)
# ============================================================

@njit
def safe_exp_arg(x):
    if x > 700.0:
        return 700.0
    if x < -700.0:
        return -700.0
    return x


@njit
def clamp01_08(x):
    if x < 0.0:
        return 0.0
    if x > 0.8:
        return 0.8
    return x


@njit
def cal_diapause(dwk, a1, b1, a2, b2, enable):
    if enable == 0:
        return 0.0, 0.0
    exp_arg0 = safe_exp_arg(-a1 * (b1 - dwk))
    exp_arg1 = safe_exp_arg(a2 * (b2 - dwk))
    z0 = 1.0 / (1.0 + math.exp(exp_arg0))
    z1 = 1.0 / (1.0 + math.exp(exp_arg1))
    return z0, z1


@njit
def cal_carrying_capacity(pn, cc, enable):
    if enable == 0:
        ccl = cc
    else:
        ccl = cc * pn
    if ccl <= 0.0:
        ccl = 1.0e-6
    return ccl


@njit
def cal_development_rate(ta, tw):
    if tw > D0_TW_MIN:
        d0 = 0.507 * math.exp(-1.0 * ((tw - 30.85) / 12.82))
        if d0 < 0.0:
            d0 = 0.0
    else:
        d0 = 1.0 / 60.0
    if tw > D1_TW_MIN:
        d1 = 0.1727 * math.exp(-1.0 * ((tw - 28.40) / 10.20))
        if d1 < 0.0:
            d1 = 0.0
    else:
        d1 = 1.0 / 60.0
    if tw > D2_TW_MIN:
        d2 = 0.602 * math.exp(-1.0 * ((tw - 34.29) / 15.07))
        if d2 < 0.0:
            d2 = 0.0
    else:
        d2 = 1.0 / 60.0
    if ta > D3_TA_MIN:
        d3 = -15.837 + 1.2897 * ta - 0.0163 * (ta * ta)
        if d3 < 0.0:
            d3 = 0.0
    else:
        d3 = 0.0
    return d0, d1, d2, d3


@njit
def cal_mortality_rate(ta, tw):
    if M0_TW_MIN < tw < 31.12:
        m0 = 0.05
    else:
        m0 = 1.0
    if M1_TW_MIN < tw < M1_TW_MAX:
        denom1 = -0.1305 * (tw * tw) + 3.868 * tw + 30.83
        denom2 = -0.1502 * (tw * tw) + 5.057 * tw + 3.517
        if abs(denom1) < 1.0e-12:
            m1 = 1.0
        else:
            m1 = 1.0 / denom1
        if abs(denom2) < 1.0e-12:
            m2 = 1.0
        else:
            m2 = 1.0 / denom2
        if m1 < 0.0:
            m1 = 0.0
        if m2 < 0.0:
            m2 = 0.0
    else:
        m1 = 1.0
        m2 = 1.0
    if M3_TA_MIN < ta < M3_TA_MAX:
        denom3 = -0.1921 * (ta * ta) + 8.147 * ta - 22.98
        if abs(denom3) < 1.0e-12:
            m3 = 1.0
        else:
            m3 = 1.0 / denom3
        if m3 < 0.0:
            m3 = 0.0
    else:
        m3 = 1.0
    if ta <= M3_TA_MIN:
        m3 = 1.0
    return m0, m1, m2, m3


@njit
def cal_washout_mortality(roff, pf, pf_scalar, enable, enable_new):
    if enable == 0:
        return 0.0
    if roff <= pf:
        return 0.0
    if enable_new == 0:
        return 1.0
    mpf = pf_scalar * (roff - pf)
    if mpf < 0.0:
        return 0.0
    if mpf > 0.8:
        return 0.8
    return mpf


@njit
def unpack_params(x):
    b1 = x[0]
    b2 = x[1]
    p_ccl = x[2]
    a1l = x[3]
    a2l = x[4]
    pf = x[5]
    p_pf_scalar_l = x[6]
    p_dis_l = x[7]
    m_mult_egg = x[8]
    m_mult_larva = x[9]
    m_mult_pupa = x[10]
    m_mult_adult = x[11]
    m_mult_diapause = x[12]
    a1 = 10.0 ** a1l
    a2 = 10.0 ** a2l
    cc = 10.0 ** p_ccl
    pf_scalar = 10.0 ** p_pf_scalar_l
    dis = 10.0 ** p_dis_l
    return (b1, b2, a1, a2, cc, pf, pf_scalar, dis,
            m_mult_egg, m_mult_larva, m_mult_pupa, m_mult_adult, m_mult_diapause)


# ============================================================
# DYNAMICS (species-specific) + STEP (integrator-specific)
# ============================================================

@njit
def dynamics_ae(y, z0, z1, ccl, d0, d1, d2, d3, m0, m1, m2, m3, m4, mpf):
    # Egg-diapause structure (Aedes albopictus).
    x0 = y[0]; x1 = y[1]; x2 = y[2]; x3 = y[3]; x4 = y[4]
    dy = np.zeros(5, dtype=np.float64)
    dy[0] = (1.0 - z0) * d3 * x3 - (m0 + d0 + mpf) * x0 + z1 * x4
    dy[1] = d0 * x0 - ((m1 + d1 + mpf) + (x1 / ccl)) * x1
    dy[2] = d1 * x1 - (m2 + d2 + mpf) * x2
    dy[3] = d2 * x2 - m3 * x3
    dy[4] = z0 * d3 * x3 - (z1 + m4) * x4
    return dy


@njit
def dynamics_cx(y, z0, z1, ccl, d0, d1, d2, d3, m0, m1, m2, m3, m4, mpf):
    # Adult-diapause structure (Culex pipiens).
    x0 = y[0]; x1 = y[1]; x2 = y[2]; x3 = y[3]; x4 = y[4]
    dy = np.zeros(5, dtype=np.float64)
    dy[0] = d3 * x3 - (m0 + d0 + mpf) * x0
    dy[1] = d0 * x0 - ((m1 + d1 + mpf) + (x1 / ccl)) * x1
    dy[2] = d1 * x1 - (m2 + d2 + mpf) * x2
    dy[3] = d2 * x2 - (m3 + z0) * x3 + z1 * x4
    dy[4] = z0 * x3 - (z1 + m4) * x4
    return dy


# Integrator steps. Each advances stage_n in place by one sub-daily DELTA_T.
# is_ae selects the diapause structure; use_rk4 selects the numerical scheme.

@njit
def _deriv(y, is_ae, z0, z1, ccl, d0, d1, d2, d3, m0, m1, m2, m3, m4, mpf):
    if is_ae:
        return dynamics_ae(y, z0, z1, ccl, d0, d1, d2, d3, m0, m1, m2, m3, m4, mpf)
    else:
        return dynamics_cx(y, z0, z1, ccl, d0, d1, d2, d3, m0, m1, m2, m3, m4, mpf)


@njit
def step_once(stage_n, is_ae, use_rk4, z0, z1, ccl, d0, d1, d2, d3, m0, m1, m2, m3, m4, mpf):
    if use_rk4:
        k1 = _deriv(stage_n, is_ae, z0, z1, ccl, d0, d1, d2, d3, m0, m1, m2, m3, m4, mpf)
        y2 = np.empty(5, dtype=np.float64)
        for i in range(5):
            y2[i] = stage_n[i] + 0.5 * DELTA_T * k1[i]
            if y2[i] < 0.0:
                y2[i] = 0.0
        k2 = _deriv(y2, is_ae, z0, z1, ccl, d0, d1, d2, d3, m0, m1, m2, m3, m4, mpf)
        y3 = np.empty(5, dtype=np.float64)
        for i in range(5):
            y3[i] = stage_n[i] + 0.5 * DELTA_T * k2[i]
            if y3[i] < 0.0:
                y3[i] = 0.0
        k3 = _deriv(y3, is_ae, z0, z1, ccl, d0, d1, d2, d3, m0, m1, m2, m3, m4, mpf)
        y4 = np.empty(5, dtype=np.float64)
        for i in range(5):
            y4[i] = stage_n[i] + DELTA_T * k3[i]
            if y4[i] < 0.0:
                y4[i] = 0.0
        k4 = _deriv(y4, is_ae, z0, z1, ccl, d0, d1, d2, d3, m0, m1, m2, m3, m4, mpf)
        for i in range(5):
            stage_n[i] = stage_n[i] + (DELTA_T / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i])
            if stage_n[i] < 0.0:
                stage_n[i] = 0.0
    else:
        # Explicit Euler.
        k1 = _deriv(stage_n, is_ae, z0, z1, ccl, d0, d1, d2, d3, m0, m1, m2, m3, m4, mpf)
        for i in range(5):
            stage_n[i] = stage_n[i] + DELTA_T * k1[i]
            if stage_n[i] < 0.0:
                stage_n[i] = 0.0


# ============================================================
# SIMULATION CORE (njit): weekly series and direct log-likelihood
# ============================================================

@njit
def _run_warmup(stage_n, is_ae, use_rk4,
                ta_wup, tw_wup, dwk_wup, smoi_wup, roff_wup,
                a1, b1, a2, b2, cc, pf, pf_scalar,
                m_mult_egg, m_mult_larva, m_mult_pupa, m_mult_adult, m_mult_diapause,
                enable_diapause, enable_ccl, enable_wash, enable_new_wash):
    for day in range(ta_wup.shape[0]):
        ta_val = ta_wup[day]
        tw_val = tw_wup[day]
        dwk_val = dwk_wup[day]
        smoi_val = smoi_wup[day] / 100.0
        roff_val = roff_wup[day]
        z0, z1 = cal_diapause(dwk_val, a1, b1, a2, b2, enable_diapause)
        ccl = cal_carrying_capacity(smoi_val, cc, enable_ccl)
        d0, d1, d2, d3 = cal_development_rate(ta_val, tw_val)
        m0, m1, m2, m3 = cal_mortality_rate(ta_val, tw_val)
        m0 = clamp01_08(m0 * m_mult_egg)
        m1 = clamp01_08(m1 * m_mult_larva)
        m2 = clamp01_08(m2 * m_mult_pupa)
        m3 = clamp01_08(m3 * m_mult_adult)
        if is_ae:
            m4 = clamp01_08(m0 * m_mult_diapause)
        else:
            m4 = clamp01_08(M4_BASE_CX * m_mult_diapause)
        mpf = cal_washout_mortality(roff_val, pf, pf_scalar, enable_wash, enable_new_wash)
        for _ in range(DAY_LENGTH):
            step_once(stage_n, is_ae, use_rk4, z0, z1, ccl, d0, d1, d2, d3, m0, m1, m2, m3, m4, mpf)


@njit
def simulate_loglik_direct(
    ta_wup, tw_wup, dwk_wup, smoi_wup, roff_wup,
    ta, tw, dwk, smoi, roff, data_year, obs,
    x, is_ae, use_rk4,
    enable_diapause, enable_ccl, enable_wash, enable_new_wash,
):
    (b1, b2, a1, a2, cc, pf, pf_scalar, dis,
     m_mult_egg, m_mult_larva, m_mult_pupa, m_mult_adult, m_mult_diapause) = unpack_params(x)

    stage_n = np.zeros(5, dtype=np.float64)
    stage_n[3] = 100.0
    stage_n[4] = 100.0

    _run_warmup(stage_n, is_ae, use_rk4,
                ta_wup, tw_wup, dwk_wup, smoi_wup, roff_wup,
                a1, b1, a2, b2, cc, pf, pf_scalar,
                m_mult_egg, m_mult_larva, m_mult_pupa, m_mult_adult, m_mult_diapause,
                enable_diapause, enable_ccl, enable_wash, enable_new_wash)

    stage_sum = np.zeros(5, dtype=np.float64)
    y = 0
    day_count_for_weekave = 1
    week_count = 0
    year_day_count = 0
    ll = 0.0

    for day in range(ta.shape[0]):
        ta_val = ta[day]
        tw_val = tw[day]
        dwk_val = dwk[day]
        smoi_val = smoi[day] / 100.0
        roff_val = roff[day]
        z0, z1 = cal_diapause(dwk_val, a1, b1, a2, b2, enable_diapause)
        ccl = cal_carrying_capacity(smoi_val, cc, enable_ccl)
        d0, d1, d2, d3 = cal_development_rate(ta_val, tw_val)
        m0, m1, m2, m3 = cal_mortality_rate(ta_val, tw_val)
        m0 = clamp01_08(m0 * m_mult_egg)
        m1 = clamp01_08(m1 * m_mult_larva)
        m2 = clamp01_08(m2 * m_mult_pupa)
        m3 = clamp01_08(m3 * m_mult_adult)
        if is_ae:
            m4 = clamp01_08(m0 * m_mult_diapause)
        else:
            m4 = clamp01_08(M4_BASE_CX * m_mult_diapause)
        mpf = cal_washout_mortality(roff_val, pf, pf_scalar, enable_wash, enable_new_wash)
        for _ in range(DAY_LENGTH):
            step_once(stage_n, is_ae, use_rk4, z0, z1, ccl, d0, d1, d2, d3, m0, m1, m2, m3, m4, mpf)
            for k in range(5):
                stage_sum[k] += stage_n[k]
        if day_count_for_weekave % 7 == 0:
            if week_count < WEEK_LENGTH:
                ev = dis * (stage_sum[3] / (DAY_LENGTH * 7.0))
                ov = obs[week_count]
                if ov != -999.0:
                    if ev < 1.0e-6:
                        ev = 1.0e-6
                    lf = 0.0
                    nn = int(ov)
                    for i in range(1, nn + 1):
                        lf += math.log(i)
                    ll += ov * math.log(ev) - ev - lf
            week_count += 1
            for k in range(5):
                stage_sum[k] = 0.0
        day_count_for_weekave += 1
        year_day_count += 1
        if y < data_year.shape[0] and year_day_count >= data_year[y]:
            y += 1
            year_day_count = 0
            day_count_for_weekave = 1
            for k in range(5):
                stage_sum[k] = 0.0
            if y >= data_year.shape[0]:
                break
        if week_count >= WEEK_LENGTH:
            break
    return ll


@njit
def simulate_est_week(
    ta_wup, tw_wup, dwk_wup, smoi_wup, roff_wup,
    ta, tw, dwk, smoi, roff, data_year,
    x, is_ae, use_rk4,
    enable_diapause, enable_ccl, enable_wash, enable_new_wash,
):
    (b1, b2, a1, a2, cc, pf, pf_scalar, dis,
     m_mult_egg, m_mult_larva, m_mult_pupa, m_mult_adult, m_mult_diapause) = unpack_params(x)

    stage_n = np.zeros(5, dtype=np.float64)
    stage_n[3] = 100.0
    stage_n[4] = 100.0

    _run_warmup(stage_n, is_ae, use_rk4,
                ta_wup, tw_wup, dwk_wup, smoi_wup, roff_wup,
                a1, b1, a2, b2, cc, pf, pf_scalar,
                m_mult_egg, m_mult_larva, m_mult_pupa, m_mult_adult, m_mult_diapause,
                enable_diapause, enable_ccl, enable_wash, enable_new_wash)

    est_week = np.zeros(WEEK_LENGTH, dtype=np.float64)
    stage_sum = np.zeros(5, dtype=np.float64)
    y = 0
    day_count_for_weekave = 1
    week_count = 0
    year_day_count = 0

    for day in range(ta.shape[0]):
        ta_val = ta[day]
        tw_val = tw[day]
        dwk_val = dwk[day]
        smoi_val = smoi[day] / 100.0
        roff_val = roff[day]
        z0, z1 = cal_diapause(dwk_val, a1, b1, a2, b2, enable_diapause)
        ccl = cal_carrying_capacity(smoi_val, cc, enable_ccl)
        d0, d1, d2, d3 = cal_development_rate(ta_val, tw_val)
        m0, m1, m2, m3 = cal_mortality_rate(ta_val, tw_val)
        m0 = clamp01_08(m0 * m_mult_egg)
        m1 = clamp01_08(m1 * m_mult_larva)
        m2 = clamp01_08(m2 * m_mult_pupa)
        m3 = clamp01_08(m3 * m_mult_adult)
        if is_ae:
            m4 = clamp01_08(m0 * m_mult_diapause)
        else:
            m4 = clamp01_08(M4_BASE_CX * m_mult_diapause)
        mpf = cal_washout_mortality(roff_val, pf, pf_scalar, enable_wash, enable_new_wash)
        for _ in range(DAY_LENGTH):
            step_once(stage_n, is_ae, use_rk4, z0, z1, ccl, d0, d1, d2, d3, m0, m1, m2, m3, m4, mpf)
            for k in range(5):
                stage_sum[k] += stage_n[k]
        if day_count_for_weekave % 7 == 0:
            if week_count < WEEK_LENGTH:
                est_week[week_count] = dis * (stage_sum[3] / (DAY_LENGTH * 7.0))
            week_count += 1
            for k in range(5):
                stage_sum[k] = 0.0
        day_count_for_weekave += 1
        year_day_count += 1
        if y < data_year.shape[0] and year_day_count >= data_year[y]:
            y += 1
            year_day_count = 0
            day_count_for_weekave = 1
            for k in range(5):
                stage_sum[k] = 0.0
            if y >= data_year.shape[0]:
                break
        if week_count >= WEEK_LENGTH:
            break
    return est_week


# ============================================================
# STATISTICS
# ============================================================

def calculate_stats(obs, est, logL, k):
    obs = np.asarray(obs, dtype=float)
    est = np.asarray(est, dtype=float)
    valid = (obs >= 0.0) & np.isfinite(obs) & np.isfinite(est)
    n = int(np.sum(valid))
    out = {
        "n": n,
        "k": int(k),
        "logL": float(logL),
        "NLL": float(-logL),
        "AIC": float(2 * k - 2 * logL),
        "BIC": float(k * math.log(n) - 2 * logL) if n > 0 else float("nan"),
    }
    if n > k + 1:
        out["AICc"] = float(out["AIC"] + (2 * k * (k + 1)) / (n - k - 1))
    else:
        out["AICc"] = float("nan")
    if n == 0:
        return out
    y = obs[valid]
    yhat = est[valid]
    ss_res = float(np.sum((y - yhat) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    out["R2_Linear"] = float(1.0 - ss_res / ss_tot) if ss_tot > 0 else float("nan")
    out["RMSE_Linear"] = float(np.sqrt(np.mean((y - yhat) ** 2)))
    ly = np.log10(y + 1.0)
    lyhat = np.log10(yhat + 1.0)
    ss_res_l = float(np.sum((ly - lyhat) ** 2))
    ss_tot_l = float(np.sum((ly - np.mean(ly)) ** 2))
    out["R2_Log10"] = float(1.0 - ss_res_l / ss_tot_l) if ss_tot_l > 0 else float("nan")
    out["RMSE_Log10"] = float(np.sqrt(np.mean((ly - lyhat) ** 2)))
    out["MAE"] = float(np.mean(np.abs(yhat - y)))
    out["Bias"] = float(np.mean(yhat - y))
    if n >= 3 and SCIPY_AVAILABLE:
        try:
            out["Pearson_r_Linear"] = float(pearsonr(y, yhat)[0])
        except Exception:
            out["Pearson_r_Linear"] = float("nan")
        try:
            out["Pearson_r_Log10"] = float(pearsonr(ly, lyhat)[0])
        except Exception:
            out["Pearson_r_Log10"] = float("nan")
        try:
            out["Spearman_r_Linear"] = float(spearmanr(y, yhat).correlation)
        except Exception:
            out["Spearman_r_Linear"] = float("nan")
        try:
            out["Spearman_r_Log10"] = float(spearmanr(ly, lyhat).correlation)
        except Exception:
            out["Spearman_r_Log10"] = float("nan")
    return out


# ============================================================
# ESTIMATION
# ============================================================

def get_active_indices(flags):
    new_wash, dis_on, mort_on = flags
    active = list(CORE_IDX)
    if new_wash:
        active += PF_SCALAR_IDX
    if dis_on:
        active += DIS_IDX
    if mort_on:
        active += MORT_IDX
    return np.array(active, dtype=np.int64)


def make_fixed_x(bounds, flags):
    new_wash, dis_on, mort_on = flags
    x = np.array([(lo + hi) / 2.0 for lo, hi in bounds], dtype=np.float64)
    if not new_wash:
        x[6] = FIXED_PF_SCALAR_L
    if not dis_on:
        x[7] = FIXED_DIS_L
    if not mort_on:
        x[8:13] = FIXED_M_MULT
    return x


class Objective:
    def __init__(self, arrays, fixed_x, active_idx, is_ae, use_rk4, enables, label):
        self.arrays = arrays
        self.fixed_x = fixed_x
        self.active_idx = active_idx
        self.is_ae = is_ae
        self.use_rk4 = use_rk4
        self.enables = enables  # (diapause, ccl, wash, new_wash)
        self.label = label
        self.count = 0
        self.best_nll = float("inf")
        self.best_x = None
        self.t0 = time.time()

    def __call__(self, theta):
        theta = np.asarray(theta, dtype=np.float64)
        self.count += 1
        if not np.all(np.isfinite(theta)):
            return MAX_PENALTY
        x = self.fixed_x.copy()
        for j, idx in enumerate(self.active_idx):
            x[idx] = theta[j]
        a = self.arrays
        ed, ec, ew, enw = self.enables
        try:
            ll = float(simulate_loglik_direct(
                a[0], a[1], a[2], a[3], a[5],
                a[6], a[7], a[8], a[9], a[11], a[13], a[12],
                x, self.is_ae, self.use_rk4, ed, ec, ew, enw,
            ))
        except Exception:
            return MAX_PENALTY
        if not math.isfinite(ll):
            return MAX_PENALTY
        nll = -ll
        if nll < self.best_nll:
            self.best_nll = nll
            self.best_x = x.copy()
        if self.count % PRINT_EVERY_OBJECTIVE == 0:
            print(f"[{self.label}] eval={self.count}, bestNLL={self.best_nll:.3f}, "
                  f"elapsed={time.time()-self.t0:.1f}s", flush=True)
        return nll


def run_one(species_key, integrator, config_key, arrays,
            de_maxiter, de_popsize):
    info = SPECIES_INFO[species_key]
    is_ae = (species_key == "Ae")
    use_rk4 = (integrator == "R")
    flags = CONFIG_FLAGS[config_key]
    bounds_full = info["bounds"]

    enable_diapause = 1
    enable_ccl = 1
    enable_wash = 1
    enable_new_wash = 1 if flags[0] else 0
    enables = (enable_diapause, enable_ccl, enable_wash, enable_new_wash)

    active_idx = get_active_indices(flags)
    fixed_x = make_fixed_x(bounds_full, flags)
    active_bounds = [(float(bounds_full[i, 0]), float(bounds_full[i, 1])) for i in active_idx]
    k = len(active_idx)

    label = build_label(config_key, integrator)
    full_name = build_full_name(config_key, integrator, info["diapause"])

    print("\n" + "=" * 72)
    print(f"Running: {label}  ({info['pretty']})")
    print(f"  {full_name}")
    print(f"  active parameters (k={k}): "
          + ", ".join(ALL_PARAM_NAMES[i] for i in active_idx))
    print("=" * 72, flush=True)

    if not SCIPY_AVAILABLE:
        raise ImportError("scipy is required for Differential Evolution.")

    # JIT warm-up / compile.
    theta0 = np.array([(lo + hi) / 2.0 for lo, hi in active_bounds], dtype=np.float64)
    x0 = fixed_x.copy()
    for j, idx in enumerate(active_idx):
        x0[idx] = theta0[j]
    a = arrays
    _ = simulate_loglik_direct(a[0], a[1], a[2], a[3], a[5],
                               a[6], a[7], a[8], a[9], a[11], a[13], a[12],
                               x0, is_ae, use_rk4, *enables)
    _ = simulate_est_week(a[0], a[1], a[2], a[3], a[5],
                          a[6], a[7], a[8], a[9], a[11], a[13],
                          x0, is_ae, use_rk4, *enables)

    obj = Objective(arrays, fixed_x, active_idx, is_ae, use_rk4, enables, label)
    t0 = time.time()
    result = differential_evolution(
        obj, bounds=active_bounds,
        maxiter=int(de_maxiter), popsize=int(de_popsize),
        tol=DE_TOL, mutation=DE_MUTATION, recombination=DE_RECOMBINATION,
        seed=DE_SEED, polish=DE_POLISH, workers=DE_WORKERS,
        updating=DE_UPDATING, init="latinhypercube",
    )
    elapsed = time.time() - t0

    # Best parameter vector (full length 13).
    xbest = fixed_x.copy()
    for j, idx in enumerate(active_idx):
        xbest[idx] = result.x[j]
    if obj.best_x is not None and obj.best_nll <= -(-result.fun):
        # prefer tracked best if strictly better
        if obj.best_nll < result.fun:
            xbest = obj.best_x.copy()
    best_nll = min(result.fun, obj.best_nll)
    best_logL = -best_nll

    est = simulate_est_week(a[0], a[1], a[2], a[3], a[5],
                            a[6], a[7], a[8], a[9], a[11], a[13],
                            xbest, is_ae, use_rk4, *enables)
    obs = arrays[12]
    stats = calculate_stats(obs, est, best_logL, k)

    return {
        "species_key": species_key,
        "integrator": integrator,
        "config_key": config_key,
        "label": label,
        "full_name": full_name,
        "flags": {"new_washout": flags[0], "s_obs_estimated": flags[1],
                  "mortality_estimated": flags[2]},
        "k": k,
        "active_parameters": [ALL_PARAM_NAMES[i] for i in active_idx],
        "xbest": xbest,
        "est": est,
        "obs": obs,
        "best_logL": best_logL,
        "best_nll": best_nll,
        "stats": stats,
        "elapsed": elapsed,
        "n_eval": obj.count,
        "scipy_success": bool(result.success),
        "scipy_message": str(result.message),
        "de_maxiter": int(de_maxiter),
        "de_popsize": int(de_popsize),
    }


# ============================================================
# OUTPUT
# ============================================================

def res_flags_tuple(res):
    f = res["flags"]
    return (f["new_washout"], f["s_obs_estimated"], f["mortality_estimated"])


def write_outputs(res):
    info = SPECIES_INFO[res["species_key"]]
    integ = "Euler" if res["integrator"] == "E" else "RK4"
    config_code = res["label"].replace("[", "_").replace("]", "").replace("-", "_")

    script_dir = os.path.dirname(os.path.abspath(__file__))
    json_name = f"{info['code']}_{res['label']}_{POINT}_params.json".replace("[", "_").replace("]", "")
    json_path = os.path.join(script_dir, json_name)

    out_dir = f"./output_paraest/{POINT}/{info['full']}/{config_code}"
    os.makedirs(out_dir, exist_ok=True)
    csv_path = os.path.join(out_dir, f"{info['code']}_{config_code}_est_week_{POINT}.csv")
    np.savetxt(csv_path, res["est"], delimiter=",")

    xbest = res["xbest"]
    active_set = set(int(i) for i in get_active_indices(res_flags_tuple(res)))

    # Build a tidy parameter table (raw value, converted value where applicable,
    # and whether the parameter was estimated or held fixed for this config).
    converted_map = {2: ("cc", 10.0 ** xbest[2]),
                     3: ("a1", 10.0 ** xbest[3]),
                     4: ("a2", 10.0 ** xbest[4]),
                     6: ("pf_scalar", 10.0 ** xbest[6]),
                     7: ("s_obs", 10.0 ** xbest[7])}
    param_rows = []
    for i, name in enumerate(ALL_PARAM_NAMES):
        conv_name, conv_val = converted_map.get(i, ("", ""))
        param_rows.append({
            "index": i,
            "parameter": name,
            "raw_value": float(xbest[i]),
            "converted_name": conv_name,
            "converted_value": (float(conv_val) if conv_name else ""),
            "status": ("estimated" if i in active_set else "fixed"),
        })
    param_df = pd.DataFrame(param_rows)
    param_csv_path = os.path.join(out_dir, f"{info['code']}_{config_code}_params_{POINT}.csv")
    param_df.to_csv(param_csv_path, index=False)
    out = {
        "label": res["label"],
        "full_name": res["full_name"],
        "species": info["pretty"],
        "species_code": info["code"],
        "diapause": info["diapause"],
        "integration": integ,
        "optimizer": "Differential Evolution",
        "config_key": res["config_key"],
        "flags": res["flags"],
        "day_length": DAY_LENGTH,
        "delta_t": DELTA_T,
        "DE_MAXITER": res["de_maxiter"],
        "DE_POPSIZE": res["de_popsize"],
        "DE_SEED": DE_SEED,
        "n_objective_evaluations": int(res["n_eval"]),
        "elapsed_sec": float(res["elapsed"]),
        "best_log_likelihood": float(res["best_logL"]),
        "best_negative_log_likelihood": float(res["best_nll"]),
        "k_estimated_parameters": res["k"],
        "active_parameters": res["active_parameters"],
        "parameters": {name: float(xbest[i]) for i, name in enumerate(ALL_PARAM_NAMES)},
        "converted_parameters": {
            "a1": float(10.0 ** xbest[3]),
            "a2": float(10.0 ** xbest[4]),
            "cc": float(10.0 ** xbest[2]),
            "pf_scalar": float(10.0 ** xbest[6]),
            "s_obs": float(10.0 ** xbest[7]),
            "dis": float(10.0 ** xbest[7]),
        },
        "stats": res["stats"],
        "scipy_result_success": res["scipy_success"],
        "scipy_result_message": res["scipy_message"],
        "csv_path": csv_path,
        "param_csv_path": param_csv_path,
        "source_design": "unified PCMP estimation program (species/integrator/config selectable, DE optimizer)",
        "generated_at": datetime.now().isoformat(timespec="seconds"),
    }
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)

    print("\nRESULT")
    print(f"  Label      : {res['label']}")
    print(f"  Full name  : {res['full_name']}")
    print(f"  logL       : {res['best_logL']:.4f}")
    s = res["stats"]
    print(f"  AIC        : {s.get('AIC', float('nan')):.3f}")
    print(f"  R2_Log10   : {s.get('R2_Log10', float('nan')):.4f}")
    print(f"  R2_Linear  : {s.get('R2_Linear', float('nan')):.4f}")
    print(f"  RMSE_Log10 : {s.get('RMSE_Log10', float('nan')):.4f}")

    print("\n  Estimated parameters (raw scale):")
    print(f"    {'parameter':20s} {'value':>14s}  {'status':9s} {'converted':>22s}")
    for row in param_rows:
        conv = ""
        if row["converted_name"]:
            conv = f"{row['converted_name']} = {row['converted_value']:.6g}"
        print(f"    {row['parameter']:20s} {row['raw_value']:>14.6f}  "
              f"{row['status']:9s} {conv:>22s}")

    print("\n  Converted parameters (natural scale):")
    print(f"    cc        = {10.0 ** xbest[2]:.6g}")
    print(f"    a1        = {10.0 ** xbest[3]:.6g}")
    print(f"    a2        = {10.0 ** xbest[4]:.6g}")
    print(f"    pf        = {xbest[5]:.6g}")
    print(f"    pf_scalar = {10.0 ** xbest[6]:.6g}")
    print(f"    s_obs     = {10.0 ** xbest[7]:.6g}   (= dis)")

    print(f"\n  Saved JSON       : {json_path}")
    print(f"  Saved weekly CSV : {csv_path}")
    print(f"  Saved param CSV  : {param_csv_path}")
    return json_path, csv_path, out


def write_summary(all_outs, species_key, integrator):
    if not all_outs:
        return
    info = SPECIES_INFO[species_key]
    script_dir = os.path.dirname(os.path.abspath(__file__))
    rows = []
    for o in all_outs:
        row = {
            "label": o["label"],
            "config_key": o["config_key"],
            "full_name": o["full_name"],
            "new_washout": o["flags"]["new_washout"],
            "s_obs_estimated": o["flags"]["s_obs_estimated"],
            "mortality_estimated": o["flags"]["mortality_estimated"],
            "k": o["k_estimated_parameters"],
            "logL": o["best_log_likelihood"],
        }
        row.update(o["stats"])
        rows.append(row)
    df = pd.DataFrame(rows)
    if "AIC" in df.columns and len(df) > 0:
        best_aic = df["AIC"].min()
        df["delta_AIC_vs_best"] = df["AIC"] - best_aic
        if (df["config_key"] == "MA").any():
            full_aic = df.loc[df["config_key"] == "MA", "AIC"].iloc[0]
            df["delta_AIC_vs_full"] = df["AIC"] - full_aic
        d = df["delta_AIC_vs_best"].to_numpy(dtype=float)
        w = np.exp(-0.5 * d)
        df["AIC_weight"] = w / np.sum(w)
        df = df.sort_values("AIC").reset_index(drop=True)
    out_csv = os.path.join(script_dir, f"{info['code']}_M{integrator}_summary.csv")
    df.to_csv(out_csv, index=False)
    print(f"\nSaved summary: {out_csv}")
    with pd.option_context("display.max_columns", None, "display.width", 200):
        cols = [c for c in ["label", "k", "logL", "AIC", "delta_AIC_vs_best",
                            "delta_AIC_vs_full", "AIC_weight", "R2_Log10",
                            "R2_Linear", "RMSE_Log10"] if c in df.columns]
        print(df[cols].to_string(index=False))


# ============================================================
# CLI / INTERACTIVE
# ============================================================

def prompt_choice(prompt, options):
    options_l = [o.lower() for o in options]
    while True:
        ans = input(f"{prompt} [{'/'.join(options)}]: ").strip()
        if ans.lower() in options_l:
            return options[options_l.index(ans.lower())]
        print("  invalid choice, try again.")


def main():
    parser = argparse.ArgumentParser(
        description="Unified PCMP estimation (species/integrator/config selectable, DE optimizer)."
    )
    parser.add_argument("--species", choices=["Ae", "Cx"], help="Ae or Cx")
    parser.add_argument("--integrator", choices=["E", "R"], help="E=Euler, R=RK4")
    parser.add_argument("--config", help="single config key (B_PCMP, MW, MS, MM, MA, "
                                         "MA-washout, MA-s_obs, MA-mortality), where MA is the "
                                         "full all-components model; OR the run mode 'all', which is "
                                         "NOT a model but a batch command that estimates every "
                                         "config in turn and writes a summary CSV (MA is included once).")
    parser.add_argument("--de-maxiter", type=int, default=DE_MAXITER)
    parser.add_argument("--de-popsize", type=int, default=DE_POPSIZE)
    args = parser.parse_args()

    species_key = args.species or prompt_choice("Species", ["Ae", "Cx"])
    integrator = args.integrator or prompt_choice("Integrator", ["E", "R"])
    if args.config:
        config_choice = args.config
    else:
        print("Single-model configs:", ", ".join(CONFIG_ORDER))
        print("  (MA = full all-components model)")
        print("Batch mode: 'all'  -> estimate every config above in turn + summary CSV")
        print("            (this is a run mode, not a model; MA is included once)")
        config_choice = input("Config or 'all': ").strip()

    info = SPECIES_INFO[species_key]
    print(f"\nLoading input data for {info['pretty']} at {POINT} ...")
    arrays = load_input_arrays(info["full"], info["obs_candidates"])
    obs = arrays[12]
    if np.sum(obs >= 0.0) == 0:
        raise RuntimeError("No valid observation data found.")
    print(f"Valid observation weeks: {int(np.sum(obs >= 0.0))}")

    if config_choice.lower() == "all":
        all_outs = []
        for ck in CONFIG_ORDER:
            res = run_one(species_key, integrator, ck, arrays,
                          args.de_maxiter, args.de_popsize)
            _, _, out = write_outputs(res)
            all_outs.append(out)
        write_summary(all_outs, species_key, integrator)
    else:
        if config_choice not in CONFIG_FLAGS:
            raise SystemExit(f"Unknown config '{config_choice}'. "
                             f"Choose from: {', '.join(CONFIG_ORDER)} or 'all'.")
        res = run_one(species_key, integrator, config_choice, arrays,
                      args.de_maxiter, args.de_popsize)
        write_outputs(res)

    print("\nDone.")


if __name__ == "__main__":
    main()