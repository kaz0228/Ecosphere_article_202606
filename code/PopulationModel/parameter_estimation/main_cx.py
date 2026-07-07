#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
main_cx.py  --  PCMP estimation with a Culex-specific thermal response.

This is the Culex-focused sibling of the unified PCMP program. It keeps the
ENTIRE estimation framework intact (both species, both integrators, the
PCMP / MM / ... configuration series, and the Differential Evolution
optimizer) and ADDS one thing: an option to re-estimate the *temperature
response* of Culex pipiens instead of inheriting the Aedes-derived
(Jia et al. 2016) thermal function unchanged.

------------------------------------------------------------------------------
Why this exists (the biological argument)
------------------------------------------------------------------------------
In the base model the temperature-dependent development and mortality functions
(`cal_development_rate`, `cal_mortality_rate`) are SHARED between species: their
shape and coefficients come from Aedes laboratory data (Jia et al. 2016). The
two species differ only in diapause structure and in the estimated mortality
multipliers. Because Ae. albopictus is subtropical (more heat-tolerant) while
Cx. pipiens is temperate (more cold-tolerant, less heat-tolerant), forcing the
Aedes thermal curve onto Culex tends to keep Culex alive through the summer and
carry too many adults into autumn -- exactly the summer/autumn mis-fit seen in
the Kumagaya panels.

No published study gives a Jia-style closed-form thermal function specifically
for Cx. pipiens. There ARE, however, experimental reports of qualitative facts:
near-complete adult mortality around ~35C, a survival optimum roughly in the
20-28C band, and an upper thermal limit lower than that of Aedes. The strategy
here is therefore the standard "borrow the functional FORM, re-estimate the
POSITION/WIDTH, then validate against the independent experimental literature":

    * The functional FORM of the Jia mortality curves is retained (each stage
      mortality is the inverse of a downward quadratic = a thermal-performance
      curve), but each curve is re-expressed about its vertex and re-estimated.
    * For the larva, pupa and adult stages (the three whose mortality comes from
      Jia/Aedes, ref [24] in Fukui et al. 2022), three interpretable knobs are
      ESTIMATED per stage from the Culex observation series:
          dTopt : shift (deg C) of the thermal optimum away from the Aedes value
                  (negative => colder optimum, expected for temperate Culex).
          fwarm : multiplier on the WARM-side curvature. >1 => mortality climbs
                  faster as it gets hot (heat-sensitive: the Culex trait).
          fcool : multiplier on the COOL-side curvature. <1 => better cold
                  tolerance (cold-tolerant: the Culex trait).
      The optimum depth m_min is held at the Aedes value (it is confounded with
      the stage mortality multiplier m_mult_<stage>, the separately-modelled
      local-climate/landscape death layer), so the overall scale stays with
      m_mult by construction and the thermal knobs only reshape the curve.
    * After fitting, the per-stage optimum and warm-side lethal temperature are
      checked against published Culex anchors (optimum ~20-28C, near-complete
      mortality ~35C). The program prints this validation block automatically.

These nine knobs are ESTIMATED ONLY FOR Culex and ONLY when thermal fitting is
switched on (config suffix "+T", or --temp-fit). With thermal fitting OFF, every
configuration -- including B_PCMP and MM -- reproduces the original base-model
behaviour bit for bit, so PCMP and MM still run exactly as before. Both Euler
and RK4 integrators are preserved.

Estimation: in addition to the existing Differential Evolution optimizer, this
program adds an MCMC (emcee) mode (--mcmc) that returns the full posterior,
convergence diagnostics, and credible intervals -- the "no-complaints" route
for the thermal re-estimation. A comparison plot of the Aedes vs Culex stage
mortality curves (larva / pupa / adult) is written automatically.

------------------------------------------------------------------------------
What this program is (unified base, retained verbatim)
------------------------------------------------------------------------------
Unified Physiology-based Climate-driven Mosquito Population (PCMP) estimation
program for two species and a hierarchy of model configurations, from the
baseline PCMP up to the full local-process model, under either time-integration
scheme, all estimated with a common Differential Evolution (DE) optimizer.

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

   Culex thermal re-estimation (this program's addition):
       Append "+T" to any config key, or pass --temp-fit, to additionally
       estimate the 9 Culex vertex thermal knobs (larva/pupa/adult x dTopt,fwarm,fcool). This only has
       an effect for --species Cx; for Ae the flag is ignored and the run is
       identical to the base model. Examples:
           --species Cx --integrator R --config MM+T
           --species Cx --integrator E --config B_PCMP+T
           --species Cx --integrator R --config MM --temp-fit   (same as MM+T)
       The batch mode 'all+T' (or --config all --temp-fit) runs every config
       with the Culex thermal parameters estimated.

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

# ------------------------------------------------------------
# SITE CONFIGURATION  (ported from main_saitama.py)
# ------------------------------------------------------------
# Each site specifies where its climate and observation data live and the
# estimation / warm-up year ranges. The climate point may differ from the
# observation point (Saitama observations are paired with Kumagaya climate).
# Switch with --site on the command line; the model itself is identical across
# sites, so --site tokyo reproduces the Tokyo results exactly.
SITE_CONFIG = {
    "tokyo": {
        "obs_point": "tokyo",        # measurement_data/<obs_point>/...
        "climate_point": "tokyo",    # input_data/<climate_point>/...
        "start_year": 2003,
        "year_length": 11,
        "warmup_start_year": 1991,
        "warmup_end_year": 2002,
    },
    "saitama": {
        "obs_point": "saitama",
        "climate_point": "kumagaya",     # Saitama observations use Kumagaya climate
        "start_year": 2014,
        "year_length": 6,                # 2014-2019 inclusive
        "warmup_start_year": 2002,
        "warmup_end_year": 2013,
    },
}

# Active site globals. Default to tokyo; overridden by apply_site_config() once
# --site is parsed. They stay module-level so the njit kernels that close over
# WEEK_LENGTH keep working. apply_site_config() MUST be called before the njit
# kernels are first compiled (i.e. early in main()), because numba bakes the
# value of WEEK_LENGTH into the compiled kernel.
POINT = "tokyo"
OBSPOINT = "tokyo"
CLIMATE_POINT = "tokyo"

START_YEAR = 2003
YEAR_LENGTH = 11
WEEK_LENGTH = 52 * YEAR_LENGTH
WARMUP_START_YEAR = 1991
WARMUP_END_YEAR = 2002

# When True, the photoperiod parameters (b1, b2, a1l, a2l) are held at the
# Tokyo-estimated values instead of being re-estimated. Set via --fix-photoperiod;
# only meaningful at Saitama (transfer the Tokyo diapause timing, re-estimate
# only landscape-dependent parameters and, optionally, the Cx thermal curves).
FIX_PHOTOPERIOD = False
FIX_PHOTOPERIOD_SPECIES = None   # "Ae" or "Cx"; set in run_one / run_one_mcmc


def apply_site_config(site):
    """Set the module-level site globals from SITE_CONFIG[site]. Call this BEFORE
    any njit kernel is compiled so WEEK_LENGTH is baked in correctly."""
    global POINT, OBSPOINT, CLIMATE_POINT, START_YEAR, YEAR_LENGTH
    global WEEK_LENGTH, WARMUP_START_YEAR, WARMUP_END_YEAR
    if site not in SITE_CONFIG:
        raise ValueError(f"Unknown site '{site}'. Choices: {list(SITE_CONFIG)}")
    cfg = SITE_CONFIG[site]
    OBSPOINT = cfg["obs_point"]
    POINT = cfg["obs_point"]          # kept for backward-compat references
    CLIMATE_POINT = cfg["climate_point"]
    START_YEAR = cfg["start_year"]
    YEAR_LENGTH = cfg["year_length"]
    WEEK_LENGTH = 52 * YEAR_LENGTH
    WARMUP_START_YEAR = cfg["warmup_start_year"]
    WARMUP_END_YEAR = cfg["warmup_end_year"]

# Sub-daily time step (shared by Euler and RK4 so that the two schemes are
# compared at the same temporal discretization).
DAY_LENGTH = 50
DELTA_T = 1.0 / DAY_LENGTH

# Differential Evolution settings (production defaults).
# These are LIGHT exploratory defaults chosen for speed. For final/production
# estimates, raise --de-maxiter (e.g. 600-800) and consider adding a seed.
# maxiter*popsize*D function evaluations: at maxiter=250, popsize=12, D~15 this
# is ~45k per seed, which runs quickly while still resolving most configs.
DE_MAXITER = 250
DE_POPSIZE = 12
DE_TOL = 1.0e-6
DE_SEEDS = (123,)           # single seed for speed; add seeds for final runs
DE_POLISH = False           # no local polish (faster; matches single-ML-trajectory aim)
DE_WORKERS =  -1
DE_UPDATING = "deferred"
DE_MUTATION = (0.5, 1.0)
DE_RECOMBINATION = 0.7

PRINT_EVERY_OBJECTIVE = 500
MAX_PENALTY = 1.0e12

# Fixed values used when a component is switched off.
# (log10 scale where the parameter is estimated on log10 scale.)
FIXED_DIS_L = -2.0          # s_obs = 10**(-2) = 0.01  (original PCMP fixed value)
FIXED_M_MULT_L = 0.0        # mortality multipliers neutral: 10**0 = 1.0
FIXED_PF_SCALAR_L = -2.0    # proportional-washout intensity, only used by linear type
FIXED_PFMAX_L = 0.0         # Hill washout max instantaneous rate: 10**0 = 1.0 /day
FIXED_HILL_L = math.log10(2.0)  # Hill exponent default h = 2 (matches manuscript P^2 form)


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


# --------------------------------------------------------------------------
# Culex thermal-response re-parameterisation (this program's core addition).
#
# Background (sources confirmed from Fukui et al. 2022, PLOS ONE, Table 1):
#   * The stage mortality functions m_L, m_P, m_A all come from reference [24]
#     = Jia et al. 2016 (a climate-driven mechanistic model of *Aedes
#     albopictus*). Each is the inverse of a downward quadratic in temperature:
#         m_stage(T) = 1 / (a*T^2 + b*T + c),  inside a thermal window.
#   * The egg mortality m_E comes from a DIFFERENT source [35] (Ewing 2016,
#     temperate mosquitoes), so it is NOT Aedes-derived and is left untouched.
#
# Because m_L, m_P, m_A are Aedes physiology, applying them unchanged to the
# temperate, heat-sensitive *Culex pipiens* is the suspected cause of the
# summer/autumn mis-fit. We therefore RE-ESTIMATE these three curves for Culex.
#
# Re-parameterisation (vertex form, ecologically interpretable):
#   The Jia quadratic 1/(aT^2+bT+c) is rewritten about its vertex as a U-shaped
#   mortality curve:
#         m(T) = 1 / ( 1/m_min  -  k(T) * (T - T_opt)^2 )
#   with an ASYMMETRIC curvature:
#         k(T) = k_cool   for T < T_opt   (cold side)
#         k(T) = k_warm   for T > T_opt   (warm side)
#   Interpretable parameters per stage:
#         T_opt   : temperature of MINIMUM mortality (thermal optimum)
#         k_warm  : steepness of the high-temperature rise. LARGER k_warm =>
#                   dies faster as it gets hot  (heat-sensitive; the Culex trait)
#         k_cool  : steepness of the low-temperature rise. SMALLER k_cool =>
#                   tolerates cold better        (cold-tolerant; the Culex trait)
#   m_min (the depth of the survival optimum) is HELD at the Jia value and NOT
#   estimated, because it is mathematically confounded with the stage mortality
#   multiplier m_mult_<stage> (the local-climate / landscape death layer that is
#   deliberately kept separate). Holding m_min fixed breaks that confound by
#   construction; the overall scale is carried by m_mult as intended.
#
# For Aedes (or when thermal fitting is OFF) the vertex parameters are set to
# the values that REPRODUCE the original Jia quadratic exactly (T_opt, k from
# the published a,b,c), so the base model is recovered bit-for-bit.
#
# Jia (Aedes) coefficients a,b,c and the derived vertex constants, per stage.
# Larva m_L and Pupa m_P are water-temperature driven; Adult m_A is air-driven.
# a,b,c are the EXACT coefficients used in cal_mortality_rate below.
#   larva: -0.1305, 3.868, 30.83  -> T_opt=14.82, m_min=0.01680, k=+0.1305
#   pupa : -0.1502, 5.057, 3.517  -> T_opt=16.83, m_min=0.02174, k=+0.1502
#   adult: -0.1921, 8.147, -22.98 -> T_opt=21.21, m_min=0.01579, k=+0.1921
# (k_warm = k_cool = -a reproduces the symmetric Jia parabola.)
AE_VERTEX = {
    # stage : (T_opt, m_min, k_symmetric)
    "larva": (14.8238, 0.016805, 0.1305),
    "pupa":  (16.8322, 0.021740, 0.1502),
    "adult": (21.2103, 0.015793, 0.1921),
}

# Flat scalar copies (njit kernels cannot index a Python dict).
AE_LARVA_TOPT = 14.8238
AE_LARVA_MMIN = 0.016805
AE_LARVA_K    = 0.1305
AE_PUPA_TOPT  = 16.8322
AE_PUPA_MMIN  = 0.021740
AE_PUPA_K     = 0.1502
AE_ADULT_TOPT = 21.2103
AE_ADULT_MMIN = 0.015793
AE_ADULT_K    = 0.1921

# Neutral Cx vertex parameters = the Aedes values above. With these the model is
# identical to the original; thermal fitting moves them away from these.
DTOPT_NEUTRAL = 0.0      # additive shift of T_opt away from the Aedes value
FWARM_NEUTRAL = 1.0      # multiplier on k_warm (1.0 = Aedes curvature)
FCOOL_NEUTRAL = 1.0      # multiplier on k_cool (1.0 = Aedes curvature)

# Hard mortality ceiling (per day) for the vertex form, matching clamp01_08.
CX_MORT_CAP = 0.8

# Literature anchors for Culex (used by the weak physiological prior / penalty
# and by the post-hoc validation block). Sources: survival optimum ~20-28 C and
# near-complete mortality ~35 C across Cx. pipiens studies.
CX_TOPT_LIT_LO = 20.0
CX_TOPT_LIT_HI = 28.0
CX_LETHAL_LIT  = 35.0
# Hard admissible bounds for the vertex optimum (penalised outside): keep the
# optimum physiological and the curve well-posed inside the thermal window.
CX_TOPT_HARD_LO = 12.0
CX_TOPT_HARD_HI = 32.0


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
    "p_m_mult_egg_l",
    "p_m_mult_larva_l",
    "p_m_mult_pupa_l",
    "p_m_mult_adult_l",
    "p_m_mult_diapause_l",
    "p_pf_max_l",     # NEW: Hill washout maximum instantaneous mortality rate (log10)
    "p_hill_l",       # NEW: Hill exponent h (log10), default h=2
    # NEW(Cx): vertex-form thermal re-parameterisation, 3 stages x 3 knobs.
    # Each (dTopt, fwarm, fcool) modifies that stage's Jia curve:
    #   dTopt  : additive shift of T_opt   (deg C; 0 = Aedes value)
    #   fwarm  : multiplier on warm-side curvature k_warm (1 = Aedes; >1 hotter-death)
    #   fcool  : multiplier on cool-side curvature k_cool (1 = Aedes; <1 cold-tolerant)
    "p_larva_dTopt", "p_larva_fwarm", "p_larva_fcool",
    "p_pupa_dTopt",  "p_pupa_fwarm",  "p_pupa_fcool",
    "p_adult_dTopt", "p_adult_fwarm", "p_adult_fcool",
]

# --------------------------------------------------------------------------
# Parameter bounds (identical for both species; see manuscript Methods).
# Rationale for each range is given inline. Bounds are deliberately WIDE but
# bounded away from biologically meaningless / numerically degenerate regions.
#
# index : parameter            scale    range            biological rationale
#   0,1 : b1, b2 (crit. photop) hours    [10.0, 15.0]     Kanto daylength spans
#                                                         ~9.8-14.6 h; critical
#                                                         photoperiod must lie in
#                                                         (or just beyond) the
#                                                         realised band to be
#                                                         identifiable.
#   2   : ccl (carrying cap.)   log10    [0.0, 5.0]       cc = 1 .. 1e5 individuals.
#   3,4 : a1l, a2l (steepness)  log10    [-1.0, 1.5]      a = 0.1 .. ~31.6 per hour;
#                                                         a>~30 is a step response,
#                                                         unidentifiable beyond.
#   5   : pf (washout thresh./  mm/day   [0.0, 300.0]     half-saturation runoff;
#         half-sat of Hill)                               must reach the suburban
#                                                         high-buffering regime
#                                                         (>140 mm) and the urban
#                                                         low threshold (~40 mm).
#   6   : pf_scalar (linear     log10    [-4.0, 0.0]      1e-4 .. 1.0 per mm; lower
#         washout intensity)                              than 1e-4 == no washout
#                                                         (degenerate with off).
#   7   : dis = s_obs           log10    [-3.0, 1.0]      1e-3 .. 10 observation
#                                                         scaling; site-specific.
#  8-12 : m_mult_{stage}        log10    [-2.0, 1.0]      0.01 .. 10 x baseline
#                                                         mortality; lower bound
#                                                         reaches the very stable
#                                                         suburban microhabitats.
#  13   : pf_max (Hill max rate) log10   [-1.0, 0.5]      0.1 .. ~3.16 per day max
#                                                         instantaneous washout
#                                                         mortality; never total
#                                                         loss (no "all killed").
#  14   : hill exponent h        log10   [0.0, 1.0]       h = 1 .. 10; h=2 default
#                                                         (matches P^2 Hill form).
#  15-23: Cx vertex thermal knobs, 3 stages x (dTopt, fwarm, fcool):
#         dTopt   degC   [-6.0, 4.0]   shift of T_opt from Aedes value. Negative
#                                      moves the optimum colder (temperate Cx);
#                                      range brackets the Ae->Cx optimum gap.
#         fwarm   ratio  [0.5, 4.0]    warm-side curvature multiplier; >1 means
#                                      mortality climbs faster at high temp
#                                      (heat-sensitive Cx). 1 = Aedes.
#         fcool   ratio  [0.3, 2.0]    cool-side curvature multiplier; <1 means
#                                      better cold tolerance. 1 = Aedes.
# --------------------------------------------------------------------------
_BOUNDS_COMMON = np.array([
    [11.0, 14.5], [11.0, 14.5],      # b1 (termination, spring), b2 (induction, autumn); b1>=b2 enforced; upper 14 < solstice 14.5h
    [0.0, 5.0],                      # ccl
    [0.5, 1.5], [0.5, 1.5],          # a1l, a2l  (a = ~3..32; sharper is unidentifiable at weekly resolution)
    [-1.0, 1.6],                     # pf (log10): 0.1 .. ~40 mm of RUNOFF (runoff-driven washout)
    [-4.0, 0.0],                     # pf_scalar_l
    [-3.0, 1.0],                     # dis_l
    [-2.0, 1.0], [-2.0, 1.0], [-2.0, 1.0], [-2.0, 1.0], [-2.0, 1.0],  # m_mult (log10) x5
    [-1.0, 1.0],                     # pf_max_l: 0.1 .. 10 /day (upper raised to test pinning)
    [0.0, 1.7],                      # hill_l: h = 1 .. ~50 (upper raised to allow step-like washout / threshold-limit)
    # Cx vertex thermal knobs (3 stages x dTopt,fwarm,fcool):
    [-6.0, 4.0], [0.5, 4.0], [0.3, 2.0],   # larva  dTopt, fwarm, fcool
    [-6.0, 4.0], [0.5, 4.0], [0.3, 2.0],   # pupa   dTopt, fwarm, fcool
    [-6.0, 4.0], [0.5, 4.0], [0.3, 2.0],   # adult  dTopt, fwarm, fcool
], dtype=np.float64)

BOUNDS_AE = _BOUNDS_COMMON.copy()
BOUNDS_CX = _BOUNDS_COMMON.copy()

# Index groups (positions in ALL_PARAM_NAMES).
CORE_IDX = [0, 1, 2, 3, 4, 5]   # b1, b2, ccl, a1l, a2l, pf
PF_SCALAR_IDX = [6]             # proportional (linear) washout intensity
DIS_IDX = [7]                   # observation-scaling factor s_obs
MORT_IDX = [8, 9, 10, 11, 12]   # stage-specific mortality multipliers (log10)
HILL_IDX = [13, 14]             # Hill washout: pf_max, hill exponent
CXTEMP_IDX = [15, 16, 17, 18, 19, 20, 21, 22, 23]  # Cx vertex thermal knobs (3 stages x 3)
# Per-stage sub-groups for convenience.
CX_LARVA_IDX = [15, 16, 17]
CX_PUPA_IDX  = [18, 19, 20]
CX_ADULT_IDX = [21, 22, 23]

# Which thermal stages are ESTIMATED when thermal fitting is on.
#   "all"   -> larva + pupa + adult (all 9 knobs)
#   "adult" -> adult only (3 knobs); larva & pupa held at neutral Aedes curves.
# The adult-only scope is the recommended default: across DE/MCMC runs the adult
# thermal response is stable and literature-consistent, whereas the water-stage
# (larva/pupa) knobs are not identified by the adult-only observations and swing
# wildly between fits. Set via --temp-stages on the command line.
TEMP_SCOPE = "all"


def cx_temp_active_idx():
    """Indices of the thermal knobs that are ESTIMATED for the current scope."""
    if TEMP_SCOPE == "adult":
        return list(CX_ADULT_IDX)
    return list(CXTEMP_IDX)


# Photoperiod / diapause-timing parameters (b1, b2, a1l, a2l). When
# --fix-photoperiod is used these are held at the Tokyo-estimated values below
# instead of being re-estimated, so transferring the calibration to Saitama
# re-estimates only the landscape-dependent (and, optionally, thermal) terms.
PHOTOPERIOD_IDX = [0, 1, 3, 4]

# Tokyo-estimated photoperiod values (internal scale used by unpack_params:
# b1, b2 are linear; a1l, a2l are log10(a)).
TOKYO_PHOTOPERIOD = {
    "Ae": {0: 14.500, 1: 13.586, 3: math.log10(31.619), 4: math.log10(31.623)},
    "Cx": {0: 13.988, 1: 12.607, 3: math.log10(9.185),  4: math.log10(3.162)},
}


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
# Each configuration is described by:
#   washout_type : "off"       no washout at all (theoretical only; not in series)
#                  "threshold" original PCMP step washout (runoff>pf -> rate 1.0)
#                  "linear"    proportional: pf_scalar*(runoff-pf), capped
#                  "hill"      saturating Hill: pf_max * R^h / (R^h + pf^h)
#                              (DEFAULT for the "new washout" configs)
#   dis_on       : estimate s_obs (True) vs fixed at 0.01 (False)
#   mort_on      : estimate mortality multipliers (True) vs fixed at 1.0
#
# All washout types are driven by ROFF (surface runoff already computed by the
# coupled water-balance model from precipitation AND soil moisture), so the
# "dry-soil vs saturated-soil" dependence of a given rainfall is carried by the
# runoff input rather than by raw precipitation.
#
# The integration scheme (Euler/RK4) is orthogonal and chosen separately.

# Default washout type for the "new washout" arm of the series.
NEW_WASHOUT_TYPE = "hill"

CONFIG_FLAGS = {
    # key            (washout_type,     dis,   mortality)
    "B_PCMP":        ("threshold",      False, False),
    "MW":            (NEW_WASHOUT_TYPE, False, False),  # washout only
    "MS":            ("threshold",      True,  False),  # s_obs only
    "MM":            ("threshold",      False, True),   # mortality only
    "MA":            (NEW_WASHOUT_TYPE, True,  True),   # all on (full)
    "MA-washout":    ("threshold",      True,  True),   # full minus new washout
    "MA-s_obs":      (NEW_WASHOUT_TYPE, False, True),   # full minus s_obs
    "MA-mortality":  (NEW_WASHOUT_TYPE, True,  False),  # full minus mortality
}

CONFIG_ORDER = [
    "B_PCMP", "MW", "MS", "MM", "MA",
    "MA-washout", "MA-s_obs", "MA-mortality",
]

# Numeric codes for washout type passed into njit kernels.
WASH_CODE = {"off": 0, "threshold": 1, "linear": 2, "hill": 3}

COMPONENT_LABELS = {
    "washout_off": "no washout",
    "washout_threshold": "threshold washout (original PCMP)",
    "washout_linear": "linear proportional washout",
    "washout_hill": "Hill saturating washout",
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
    washout_type, dis_on, mort_on = CONFIG_FLAGS[config_key]
    integ = "Euler" if integrator == "E" else "RK4"
    parts = [f"{integ} integration"]
    parts.append(COMPONENT_LABELS[f"washout_{washout_type}"])
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


def load_input_arrays(species_full, obs_candidates, point=None, obs_point=None,
                      climate_point=None, species_code=None):
    """Load warm-up + main climate and weekly observations for the active site.

    The climate point may differ from the observation point (Saitama
    observations are paired with Kumagaya climate), so climate is read from
    CLIMATE_POINT and observations from OBSPOINT unless overridden by arguments.
    Observation filenames differ by site / naming convention, so several
    candidate templates are tried per year (the short species code Ae/Cx matches
    the Saitama files; the long obs_candidates tokens match Tokyo)."""
    if obs_point is None:
        obs_point = OBSPOINT
    if climate_point is None:
        climate_point = CLIMATE_POINT
    cpt = climate_point

    ta_wup, tw_wup, dwk_wup, smoi_wup, prec_wup, roff_wup = [], [], [], [], [], []

    for yr in range(WARMUP_START_YEAR, WARMUP_END_YEAR + 1):
        fname = f"./../../input_data/{cpt}/{yr}/{cpt}{yr}_climdata2.csv"
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
        fname1 = f"./../../input_data/{cpt}/{yr}/{cpt}{yr}_climdata.csv"
        fname2 = f"./../../input_data/{cpt}/{yr}/{cpt}{yr}_climdata2.csv"
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

    # ---- observations -----------------------------------------------------
    obs_dir = f"./../../measurement_data/{obs_point}/{species_full}"
    short = species_code if species_code is not None else ""
    tokens = ([short] if short else []) + list(obs_candidates)  # short code first

    def candidate_paths(tok, yr):
        return [
            os.path.join(obs_dir, f"{tok}_{yr}_{obs_point}_obs.txt"),
            os.path.join(obs_dir, f"{obs_point}_{tok}_{yr}.tsv"),
            os.path.join(obs_dir, f"{obs_point}_{tok}_female_{yr}.tsv"),
        ]

    found_any = False
    for tok in tokens:
        for p in candidate_paths(tok, START_YEAR):
            if os.path.exists(p):
                print(f"Data file found: {p}")
                found_any = True
                break
        if found_any:
            break
    if not found_any:
        tried = [pp for tok in tokens for pp in candidate_paths(tok, START_YEAR)]
        raise FileNotFoundError(
            f"No observation data found in {obs_dir} for {START_YEAR}. Tried:\n  "
            + "\n  ".join(tried)
        )

    obs_adults_week = np.full(WEEK_LENGTH, -999.0, dtype=np.float64)
    for yr in range(START_YEAR, START_YEAR + YEAR_LENGTH):
        fname = None
        for tok in tokens:
            for p in candidate_paths(tok, yr):
                if os.path.exists(p):
                    fname = p
                    break
            if fname is not None:
                break
        if fname is None:
            print(f"Warning: observation file missing for year {yr} in {obs_dir}")
            continue
        df = pd.read_csv(fname, sep="\t", header=None, names=["week", "pop"], engine="python")
        week_plus = (yr - START_YEAR) * 52
        for _, row in df.iterrows():
            wk = int(row["week"])
            pop = float(row["pop"])
            wc = wk + week_plus - 1
            if 0 <= wc < WEEK_LENGTH:
                obs_adults_week[wc] = pop

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

@njit(cache=True, fastmath=True)
def safe_exp_arg(x):
    if x > 700.0:
        return 700.0
    if x < -700.0:
        return -700.0
    return x


@njit(cache=True, fastmath=True)
def clamp01_08(x):
    if x < 0.0:
        return 0.0
    if x > 0.8:
        return 0.8
    return x


@njit(cache=True, fastmath=True)
def cal_diapause(dwk, dchan, a1, b1, a2, b2, enable):
    """Daylength-driven diapause WITHOUT spring/autumn separation (legacy form).

    Both termination (z0) and induction (z1) are computed from the instantaneous
    daylength only; `dchan` is accepted for call-signature compatibility but is
    not used. This is the form that previously achieved the best fit (logL ~
    -4000). Spring/autumn separation was tried but worsened the fit while not
    resolving b1/b2 (early-spring termination is unobserved), so it was dropped.
    Temperature effects enter through development and mortality, not here.
    """
    if enable == 0:
        return 0.0, 0.0
    exp_arg0 = safe_exp_arg(-a1 * (b1 - dwk))
    exp_arg1 = safe_exp_arg(a2 * (b2 - dwk))
    z0 = 1.0 / (1.0 + math.exp(exp_arg0))
    z1 = 1.0 / (1.0 + math.exp(exp_arg1))
    return z0, z1


@njit(cache=True, fastmath=True)
def cal_carrying_capacity(pn, cc, enable):
    if enable == 0:
        ccl = cc
    else:
        ccl = cc * pn
    if ccl <= 0.0:
        ccl = 1.0e-6
    return ccl


@njit(cache=True, fastmath=True)
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


@njit(cache=True, fastmath=True)
def _stage_mort_vertex(T, Topt, m_min, k_cool, k_warm, lo, hi, cap):
    """One stage's mortality in vertex form.

        m(T) = 1 / ( 1/m_min  -  k * (T - Topt)^2 )
    with k = k_cool for T < Topt and k = k_warm for T >= Topt, clipped to [0, cap]
    and forced to cap outside the thermal window (lo, hi). This is the Aedes Jia
    quadratic when Topt, m_min come from the published a,b,c and k_cool = k_warm =
    -a; Culex moves Topt and makes k_warm > k_cool (heat-sensitive, cold-tolerant).

    DELIBERATE DIFFERENCE FROM THE ORIGINAL CODE AT THE EXTREME HIGH-TEMPERATURE
    EDGE: in the original main.py the raw quadratic 1/(aT^2+bT+c) goes NEGATIVE
    once T passes the upper root of the denominator, and clamp01_08 then turns
    that negative value into 0.0 -- i.e. mortality drops to ZERO at very hot
    temperatures (a non-physical artefact: "too hot to die"). Here, once the
    denominator collapses we return `cap` (maximum mortality) instead, which is
    the physically correct "heat kills" behaviour and is exactly the direction
    Culex needs. The two differ only in a narrow band at the top of each window
    (water > ~34.4 C, air > ~39 C); in the Culex temperature range this edge is
    essentially never reached, so the seasonal fit is unaffected. Elsewhere the
    vertex form reproduces the Jia quadratic to numerical precision.
    """
    if T <= lo or T >= hi:
        return cap
    if T < Topt:
        k = k_cool
    else:
        k = k_warm
    dt = T - Topt
    denom = (1.0 / m_min) - k * dt * dt
    if denom <= 1.0e-12:
        return cap
    m = 1.0 / denom
    if m < 0.0:
        return 0.0
    if m > cap:
        return cap
    return m


@njit(cache=True, fastmath=True)
def cal_mortality_rate(ta, tw,
                       lv_dTopt, lv_fwarm, lv_fcool,
                       pp_dTopt, pp_fwarm, pp_fcool,
                       ad_dTopt, ad_fwarm, ad_fcool):
    """Stage mortality rates (egg/diapause m0, larva m1, pupa m2, adult m3).

    Egg/overwinter mortality m0 keeps the ORIGINAL Jia logic (its window and the
    0.05/1.0 step are unchanged: the egg curve is not Aedes-physiology in the
    Fukui 2022 sourcing and is left alone).

    Larva (m1), pupa (m2) and adult (m3) use the VERTEX form. The nine arguments
    are the per-stage knobs (dTopt, fwarm, fcool). With dTopt=0, fwarm=1, fcool=1
    for every stage the vertex form reproduces the original Jia quadratics, so the
    Aedes / thermal-fit-OFF behaviour is recovered.

    larva/pupa are driven by water temperature tw; adult by air temperature ta.
    """
    # m0: unchanged egg/overwinter step (water temperature).
    if M0_TW_MIN < tw < 31.12:
        m0 = 0.05
    else:
        m0 = 1.0

    # Vertex parameters per stage = Aedes base shifted/scaled by the knobs.
    lv_Topt = AE_LARVA_TOPT + lv_dTopt
    pp_Topt = AE_PUPA_TOPT + pp_dTopt
    ad_Topt = AE_ADULT_TOPT + ad_dTopt

    lv_kc = AE_LARVA_K * lv_fcool
    lv_kw = AE_LARVA_K * lv_fwarm
    pp_kc = AE_PUPA_K * pp_fcool
    pp_kw = AE_PUPA_K * pp_fwarm
    ad_kc = AE_ADULT_K * ad_fcool
    ad_kw = AE_ADULT_K * ad_fwarm

    m1 = _stage_mort_vertex(tw, lv_Topt, AE_LARVA_MMIN, lv_kc, lv_kw,
                            M1_TW_MIN, M1_TW_MAX, CX_MORT_CAP)
    m2 = _stage_mort_vertex(tw, pp_Topt, AE_PUPA_MMIN, pp_kc, pp_kw,
                            M2_TW_MIN, M1_TW_MAX, CX_MORT_CAP)
    m3 = _stage_mort_vertex(ta, ad_Topt, AE_ADULT_MMIN, ad_kc, ad_kw,
                            M3_TA_MIN, M3_TA_MAX, CX_MORT_CAP)
    return m0, m1, m2, m3


@njit(cache=True, fastmath=True)
def cal_washout_mortality(drive, pf, pf_scalar, pf_max, hill, wash_code):
    """Instantaneous washout mortality rate (per day).

    The DRIVING variable is chosen by the caller according to wash_code:
      * threshold (original PCMP) is driven by surface RUNOFF (ROFF), which
        already folds soil moisture into the flow that actually occurred.
      * the new linear / Hill washout is driven by raw PRECIPITATION (PREC).
        Because the upstream water-balance model assumes a soil surface and
        does not distinguish impervious (urban) cover, runoff carries no
        site contrast here; precipitation-driven washout lets the SUSCEPTIBILITY
        parameters (pf, pf_max) carry the urban/suburban difference instead.
        pf is therefore a precipitation threshold / Hill half-saturation in mm.

    wash_code:
      0 off       -> no washout (0.0)
      1 threshold -> original PCMP: rate 1.0/day once the driver exceeds pf
                     (an all-or-nothing step; retained for the baseline only)
      2 linear    -> pf_scalar * (drive - pf), zero below pf, capped at 0.8
      3 hill      -> pf_max * drive^h / (drive^h + pf^h): smooth, saturating,
                     never reaches total loss (asymptote = pf_max), and gives a
                     genuinely different rate for e.g. 50 vs 100 mm of rainfall.
    """
    if wash_code == 0:
        return 0.0
    if wash_code == 1:                       # original PCMP threshold (step)
        if drive <= pf:
            return 0.0
        return 1.0
    if wash_code == 2:                       # linear proportional
        if drive <= pf:
            return 0.0
        mpf = pf_scalar * (drive - pf)
        if mpf < 0.0:
            return 0.0
        if mpf > 0.8:
            return 0.8
        return mpf
    # wash_code == 3: Hill saturating
    if drive <= 0.0:
        return 0.0
    rh = drive ** hill
    ph = pf ** hill
    denom = rh + ph
    if denom <= 0.0:
        return 0.0
    mpf = pf_max * rh / denom
    if mpf < 0.0:
        return 0.0
    # Hard safety cap well below total loss; pf_max itself is bounded so this
    # rarely binds, but it guarantees no single step removes the whole stage.
    if mpf > 0.95:
        return 0.95
    return mpf


@njit(cache=True, fastmath=True)
def unpack_params(x):
    b1 = x[0]
    b2 = x[1]
    p_ccl = x[2]
    a1l = x[3]
    a2l = x[4]
    pf = 10.0 ** x[5]         # pf now estimated on log10 scale (mm of driver)
    p_pf_scalar_l = x[6]
    p_dis_l = x[7]
    m_mult_egg_l = x[8]
    m_mult_larva_l = x[9]
    m_mult_pupa_l = x[10]
    m_mult_adult_l = x[11]
    m_mult_diapause_l = x[12]
    p_pf_max_l = x[13]
    p_hill_l = x[14]
    lv_dTopt = x[15]; lv_fwarm = x[16]; lv_fcool = x[17]
    pp_dTopt = x[18]; pp_fwarm = x[19]; pp_fcool = x[20]
    ad_dTopt = x[21]; ad_fwarm = x[22]; ad_fcool = x[23]
    a1 = 10.0 ** a1l
    a2 = 10.0 ** a2l
    cc = 10.0 ** p_ccl
    pf_scalar = 10.0 ** p_pf_scalar_l
    dis = 10.0 ** p_dis_l
    m_mult_egg = 10.0 ** m_mult_egg_l
    m_mult_larva = 10.0 ** m_mult_larva_l
    m_mult_pupa = 10.0 ** m_mult_pupa_l
    m_mult_adult = 10.0 ** m_mult_adult_l
    m_mult_diapause = 10.0 ** m_mult_diapause_l
    pf_max = 10.0 ** p_pf_max_l
    hill = 10.0 ** p_hill_l
    return (b1, b2, a1, a2, cc, pf, pf_scalar, dis,
            m_mult_egg, m_mult_larva, m_mult_pupa, m_mult_adult, m_mult_diapause,
            pf_max, hill,
            lv_dTopt, lv_fwarm, lv_fcool,
            pp_dTopt, pp_fwarm, pp_fcool,
            ad_dTopt, ad_fwarm, ad_fcool)


# ============================================================
# DYNAMICS (species-specific) + STEP (integrator-specific)
# ============================================================

@njit(cache=True, fastmath=True)
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


@njit(cache=True, fastmath=True)
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

@njit(cache=True, fastmath=True)
def _deriv(y, is_ae, z0, z1, ccl, d0, d1, d2, d3, m0, m1, m2, m3, m4, mpf):
    if is_ae:
        return dynamics_ae(y, z0, z1, ccl, d0, d1, d2, d3, m0, m1, m2, m3, m4, mpf)
    else:
        return dynamics_cx(y, z0, z1, ccl, d0, d1, d2, d3, m0, m1, m2, m3, m4, mpf)


@njit(cache=True, fastmath=True)
def step_once(stage_n, is_ae, use_rk4, z0, z1, ccl, d0, d1, d2, d3, m0, m1, m2, m3, m4, mpf):
    if use_rk4:
        k1 = _deriv(stage_n, is_ae, z0, z1, ccl, d0, d1, d2, d3, m0, m1, m2, m3, m4, mpf)
        y2 = np.empty(5, dtype=np.float64)
        for i in range(5):
            y2[i] = stage_n[i] + 0.5 * DELTA_T * k1[i]
        k2 = _deriv(y2, is_ae, z0, z1, ccl, d0, d1, d2, d3, m0, m1, m2, m3, m4, mpf)
        y3 = np.empty(5, dtype=np.float64)
        for i in range(5):
            y3[i] = stage_n[i] + 0.5 * DELTA_T * k2[i]
        k3 = _deriv(y3, is_ae, z0, z1, ccl, d0, d1, d2, d3, m0, m1, m2, m3, m4, mpf)
        y4 = np.empty(5, dtype=np.float64)
        for i in range(5):
            y4[i] = stage_n[i] + DELTA_T * k3[i]
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

@njit(cache=True, fastmath=True)
def _run_warmup(stage_n, is_ae, use_rk4,
                ta_wup, tw_wup, dwk_wup, smoi_wup, roff_wup, prec_wup,
                a1, b1, a2, b2, cc, pf, pf_scalar, pf_max, hill,
                m_mult_egg, m_mult_larva, m_mult_pupa, m_mult_adult, m_mult_diapause,
                enable_diapause, enable_ccl, wash_code,
                lv_dTopt, lv_fwarm, lv_fcool,
                pp_dTopt, pp_fwarm, pp_fcool,
                ad_dTopt, ad_fwarm, ad_fcool):
    for day in range(ta_wup.shape[0]):
        ta_val = ta_wup[day]
        tw_val = tw_wup[day]
        dwk_val = dwk_wup[day]
        dchan = (dwk_wup[day] - dwk_wup[day - 1]) if day > 0 else 1.0
        smoi_val = smoi_wup[day] / 100.0
        # all washout types are runoff-driven (runoff folds in soil moisture;
        # precipitation-only washout was tested and fit much worse)
        drive = roff_wup[day]
        z0, z1 = cal_diapause(dwk_val, dchan, a1, b1, a2, b2, enable_diapause)
        ccl = cal_carrying_capacity(smoi_val, cc, enable_ccl)
        d0, d1, d2, d3 = cal_development_rate(ta_val, tw_val)
        m0_raw, m1, m2, m3 = cal_mortality_rate(ta_val, tw_val,
                                                lv_dTopt, lv_fwarm, lv_fcool,
                                                pp_dTopt, pp_fwarm, pp_fcool,
                                                ad_dTopt, ad_fwarm, ad_fcool)
        m0 = clamp01_08(m0_raw * m_mult_egg)
        m1 = clamp01_08(m1 * m_mult_larva)
        m2 = clamp01_08(m2 * m_mult_pupa)
        m3 = clamp01_08(m3 * m_mult_adult)
        if is_ae:
            m4 = clamp01_08(m0_raw * m_mult_diapause)
        else:
            m4 = clamp01_08(M4_BASE_CX * m_mult_diapause)
        mpf = cal_washout_mortality(drive, pf, pf_scalar, pf_max, hill, wash_code)
        for _ in range(DAY_LENGTH):
            step_once(stage_n, is_ae, use_rk4, z0, z1, ccl, d0, d1, d2, d3, m0, m1, m2, m3, m4, mpf)


@njit(cache=True, fastmath=True)
def simulate_loglik_direct(
    ta_wup, tw_wup, dwk_wup, smoi_wup, roff_wup, prec_wup,
    ta, tw, dwk, smoi, roff, prec, data_year, obs,
    x, is_ae, use_rk4,
    enable_diapause, enable_ccl, wash_code,
):
    (b1, b2, a1, a2, cc, pf, pf_scalar, dis,
     m_mult_egg, m_mult_larva, m_mult_pupa, m_mult_adult, m_mult_diapause,
     pf_max, hill,
     lv_dTopt, lv_fwarm, lv_fcool,
     pp_dTopt, pp_fwarm, pp_fcool,
     ad_dTopt, ad_fwarm, ad_fcool) = unpack_params(x)
    if enable_diapause == 1 and b1 < b2:
        return -1.0e12

    stage_n = np.zeros(5, dtype=np.float64)
    stage_n[3] = 100.0
    stage_n[4] = 100.0

    _run_warmup(stage_n, is_ae, use_rk4,
                ta_wup, tw_wup, dwk_wup, smoi_wup, roff_wup, prec_wup,
                a1, b1, a2, b2, cc, pf, pf_scalar, pf_max, hill,
                m_mult_egg, m_mult_larva, m_mult_pupa, m_mult_adult, m_mult_diapause,
                enable_diapause, enable_ccl, wash_code,
                lv_dTopt, lv_fwarm, lv_fcool,
                pp_dTopt, pp_fwarm, pp_fcool,
                ad_dTopt, ad_fwarm, ad_fcool)

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
        dchan = (dwk[day] - dwk[day - 1]) if day > 0 else 1.0
        smoi_val = smoi[day] / 100.0
        drive = roff[day]   # all washout types runoff-driven (see note in warmup)
        z0, z1 = cal_diapause(dwk_val, dchan, a1, b1, a2, b2, enable_diapause)
        ccl = cal_carrying_capacity(smoi_val, cc, enable_ccl)
        d0, d1, d2, d3 = cal_development_rate(ta_val, tw_val)
        m0_raw, m1, m2, m3 = cal_mortality_rate(ta_val, tw_val,
                                                lv_dTopt, lv_fwarm, lv_fcool,
                                                pp_dTopt, pp_fwarm, pp_fcool,
                                                ad_dTopt, ad_fwarm, ad_fcool)
        m0 = clamp01_08(m0_raw * m_mult_egg)
        m1 = clamp01_08(m1 * m_mult_larva)
        m2 = clamp01_08(m2 * m_mult_pupa)
        m3 = clamp01_08(m3 * m_mult_adult)
        if is_ae:
            m4 = clamp01_08(m0_raw * m_mult_diapause)
        else:
            m4 = clamp01_08(M4_BASE_CX * m_mult_diapause)
        mpf = cal_washout_mortality(drive, pf, pf_scalar, pf_max, hill, wash_code)
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
                    lf = math.lgamma(ov + 1.0)
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


@njit(cache=True, fastmath=True)
def simulate_est_week(
    ta_wup, tw_wup, dwk_wup, smoi_wup, roff_wup, prec_wup,
    ta, tw, dwk, smoi, roff, prec, data_year,
    x, is_ae, use_rk4,
    enable_diapause, enable_ccl, wash_code,
):
    (b1, b2, a1, a2, cc, pf, pf_scalar, dis,
     m_mult_egg, m_mult_larva, m_mult_pupa, m_mult_adult, m_mult_diapause,
     pf_max, hill,
     lv_dTopt, lv_fwarm, lv_fcool,
     pp_dTopt, pp_fwarm, pp_fcool,
     ad_dTopt, ad_fwarm, ad_fcool) = unpack_params(x)
    if enable_diapause == 1 and b1 < b2:
        return np.zeros(WEEK_LENGTH, dtype=np.float64)

    stage_n = np.zeros(5, dtype=np.float64)
    stage_n[3] = 100.0
    stage_n[4] = 100.0

    _run_warmup(stage_n, is_ae, use_rk4,
                ta_wup, tw_wup, dwk_wup, smoi_wup, roff_wup, prec_wup,
                a1, b1, a2, b2, cc, pf, pf_scalar, pf_max, hill,
                m_mult_egg, m_mult_larva, m_mult_pupa, m_mult_adult, m_mult_diapause,
                enable_diapause, enable_ccl, wash_code,
                lv_dTopt, lv_fwarm, lv_fcool,
                pp_dTopt, pp_fwarm, pp_fcool,
                ad_dTopt, ad_fwarm, ad_fcool)

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
        dchan = (dwk[day] - dwk[day - 1]) if day > 0 else 1.0
        smoi_val = smoi[day] / 100.0
        drive = roff[day]   # all washout types runoff-driven (see note in warmup)
        z0, z1 = cal_diapause(dwk_val, dchan, a1, b1, a2, b2, enable_diapause)
        ccl = cal_carrying_capacity(smoi_val, cc, enable_ccl)
        d0, d1, d2, d3 = cal_development_rate(ta_val, tw_val)
        m0_raw, m1, m2, m3 = cal_mortality_rate(ta_val, tw_val,
                                                lv_dTopt, lv_fwarm, lv_fcool,
                                                pp_dTopt, pp_fwarm, pp_fcool,
                                                ad_dTopt, ad_fwarm, ad_fcool)
        m0 = clamp01_08(m0_raw * m_mult_egg)
        m1 = clamp01_08(m1 * m_mult_larva)
        m2 = clamp01_08(m2 * m_mult_pupa)
        m3 = clamp01_08(m3 * m_mult_adult)
        if is_ae:
            m4 = clamp01_08(m0_raw * m_mult_diapause)
        else:
            m4 = clamp01_08(M4_BASE_CX * m_mult_diapause)
        mpf = cal_washout_mortality(drive, pf, pf_scalar, pf_max, hill, wash_code)
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

def get_active_indices(flags, temp_fit=False):
    """Active (estimated) parameter indices for a config.

    flags = (washout_type, dis_on, mort_on).
    CORE_IDX always includes pf (index 5), which is the threshold for the
    threshold/linear types and the half-saturation runoff for the Hill type.

    temp_fit (Cx only): additionally estimate the Culex thermal parameters
    temp_fit (Cx only): additionally estimate the 9 Culex vertex thermal knobs
    (CXTEMP_IDX). The caller is responsible for only passing
    temp_fit=True when the species is Culex.
    """
    washout_type, dis_on, mort_on = flags
    active = list(CORE_IDX)            # includes pf (idx 5)
    if washout_type == "linear":
        active += PF_SCALAR_IDX        # estimate proportional intensity
    elif washout_type == "hill":
        active += HILL_IDX             # estimate pf_max and hill exponent
    # threshold / off: no extra washout parameters (only pf, already in CORE)
    if dis_on:
        active += DIS_IDX
    if mort_on:
        active += MORT_IDX
    if temp_fit:
        active += cx_temp_active_idx()  # all 9, or adult-only (3), per TEMP_SCOPE
    # When photoperiod is held at the Tokyo values, drop b1, b2, a1l, a2l from
    # the estimated set so only landscape-dependent (and thermal) params are fit.
    if FIX_PHOTOPERIOD:
        active = [i for i in active if i not in PHOTOPERIOD_IDX]
    return np.array(sorted(active), dtype=np.int64)


def make_fixed_x(bounds, flags, temp_fit=False):
    """Full-length parameter vector with inactive components held at their
    neutralising fixed values, active components at the bound midpoint."""
    washout_type, dis_on, mort_on = flags
    x = np.array([(lo + hi) / 2.0 for lo, hi in bounds], dtype=np.float64)
    # Washout parameters not used by the chosen type are pinned to neutral.
    if washout_type != "linear":
        x[6] = FIXED_PF_SCALAR_L
    if washout_type != "hill":
        x[13] = FIXED_PFMAX_L
        x[14] = FIXED_HILL_L
    if not dis_on:
        x[7] = FIXED_DIS_L
    if not mort_on:
        x[8:13] = FIXED_M_MULT_L     # log10 scale: 0.0 -> multiplier 1.0
    # Cx thermal vertex knobs: pin to NEUTRAL (Aedes values) for every stage
    # that is NOT being estimated. With thermal fitting off, all 9 are neutral
    # (= original Jia curves). With adult-only scope, larva & pupa stay neutral
    # and only the adult knobs are estimated.
    estimated_thermal = set(cx_temp_active_idx()) if temp_fit else set()
    for j in CXTEMP_IDX:
        if j in estimated_thermal:
            continue
        if (j - 15) % 3 == 0:              # dTopt slots
            x[j] = DTOPT_NEUTRAL
        else:                               # fwarm / fcool slots
            x[j] = FWARM_NEUTRAL            # == FCOOL_NEUTRAL == 1.0
    # Hold photoperiod at the Tokyo-estimated values when requested. These
    # indices are excluded from the active set in get_active_indices, so the
    # values written here are the ones actually used in the simulation.
    if FIX_PHOTOPERIOD and FIX_PHOTOPERIOD_SPECIES is not None:
        for idx, val in TOKYO_PHOTOPERIOD[FIX_PHOTOPERIOD_SPECIES].items():
            x[idx] = val
    return x


class Objective:
    def __init__(self, arrays, fixed_x, active_idx, is_ae, use_rk4, enables, label):
        self.arrays = arrays
        self.fixed_x = fixed_x
        self.active_idx = active_idx
        self.is_ae = is_ae
        self.use_rk4 = use_rk4
        self.enables = enables  # (diapause, ccl, wash_code)
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
        ed, ec, wc = self.enables
        try:
            ll = float(simulate_loglik_direct(
                a[0], a[1], a[2], a[3], a[5], a[4],
                a[6], a[7], a[8], a[9], a[11], a[10], a[13], a[12],
                x, self.is_ae, self.use_rk4, ed, ec, wc,
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


class ProgressCallback:
    """Called by differential_evolution once per generation, in the MAIN
    process (even when workers=-1). This is the parallel-safe way to show
    'it is running' progress: the per-evaluation print inside Objective runs
    in worker processes and does not reach the console, but the callback runs
    in the main process every generation.

    It evaluates the objective at the current best vector xk to report the
    best logL so far, and records the trajectory for an optional convergence
    plot.
    """
    def __init__(self, obj, label, report_every=1):
        self.obj = obj
        self.label = label
        self.report_every = report_every
        self.gen = 0
        self.t0 = time.time()
        self.history = []  # list of (gen, best_logL, elapsed_s)

    def __call__(self, xk, convergence=0.0):
        self.gen += 1
        nll = self.obj(np.asarray(xk, dtype=np.float64))
        best_logL = -nll
        elapsed = time.time() - self.t0
        self.history.append((self.gen, best_logL, elapsed))
        if self.gen % self.report_every == 0:
            print(f"[{self.label}] gen={self.gen:4d}  best logL={best_logL:12.3f}  "
                  f"convergence={convergence:.3e}  elapsed={elapsed:6.1f}s", flush=True)
        return False  # returning True would stop the optimization


def run_one(species_key, integrator, config_key, arrays,
            de_maxiter, de_popsize, temp_fit=False):
    info = SPECIES_INFO[species_key]
    is_ae = (species_key == "Ae")
    use_rk4 = (integrator == "R")
    flags = CONFIG_FLAGS[config_key]
    bounds_full = info["bounds"]

    # Temperature fitting is a Culex-only feature; silently ignore for Aedes so
    # that an --temp-fit batch run does not change any Aedes result.
    effective_temp_fit = bool(temp_fit) and (species_key == "Cx")
    global FIX_PHOTOPERIOD_SPECIES
    FIX_PHOTOPERIOD_SPECIES = species_key

    enable_diapause = 1
    enable_ccl = 1
    wash_code = WASH_CODE[flags[0]]
    enables = (enable_diapause, enable_ccl, wash_code)

    active_idx = get_active_indices(flags, effective_temp_fit)
    fixed_x = make_fixed_x(bounds_full, flags, effective_temp_fit)
    active_bounds = [(float(bounds_full[i, 0]), float(bounds_full[i, 1])) for i in active_idx]
    k = len(active_idx)

    label = build_label(config_key, integrator)
    if effective_temp_fit:
        label = label + ("+Ta" if TEMP_SCOPE == "adult" else "+T")
    full_name = build_full_name(config_key, integrator, info["diapause"])
    if effective_temp_fit:
        full_name = full_name + "; Culex thermal curves re-estimated (vertex form, larva/pupa/adult)"

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
    _ = simulate_loglik_direct(a[0], a[1], a[2], a[3], a[5], a[4],
                               a[6], a[7], a[8], a[9], a[11], a[10], a[13], a[12],
                               x0, is_ae, use_rk4, *enables)
    _ = simulate_est_week(a[0], a[1], a[2], a[3], a[5], a[4],
                          a[6], a[7], a[8], a[9], a[11], a[10], a[13],
                          x0, is_ae, use_rk4, *enables)

    # Multi-start DE: run independently from several seeds and keep the best
    # optimum. A single seeded run can settle in a local / degenerate optimum
    # (the source of the earlier pf_scalar~6000, s_obs~1e-5 artefacts); the
    # best-over-seeds rule guards against that.
    t0 = time.time()
    best_result = None
    best_obj = None
    n_eval_total = 0
    for seed in DE_SEEDS:
        obj = Objective(arrays, fixed_x, active_idx, is_ae, use_rk4, enables, label)
        progress = ProgressCallback(obj, f"{label}|seed{seed}", report_every=25)
        result = differential_evolution(
            obj, bounds=active_bounds,
            maxiter=int(de_maxiter), popsize=int(de_popsize),
            tol=DE_TOL, mutation=DE_MUTATION, recombination=DE_RECOMBINATION,
            seed=seed, polish=DE_POLISH, workers=DE_WORKERS,
            updating=DE_UPDATING, init="latinhypercube",
            callback=progress,
        )
        n_eval_total += int(result.nfev) if hasattr(result, "nfev") else int(obj.count)
        print(f"  [seed {seed}] logL = {-float(result.fun):.4f}", flush=True)
        if (best_result is None) or (float(result.fun) < float(best_result.fun)):
            best_result = result
            best_obj = obj
    result = best_result
    obj = best_obj
    elapsed = time.time() - t0
    print(f"  best-over-seeds logL = {-float(result.fun):.4f}", flush=True)

    # Best parameter vector (full length 15).
    # NOTE: differential_evolution returns the best solution of the run in
    # result.x / result.fun; we rely on these directly (correct in both serial
    # and parallel workers=-1 modes).
    xbest = fixed_x.copy()
    for j, idx in enumerate(active_idx):
        xbest[idx] = result.x[j]
    best_nll = float(result.fun)
    best_logL = -best_nll

    est = simulate_est_week(a[0], a[1], a[2], a[3], a[5], a[4],
                            a[6], a[7], a[8], a[9], a[11], a[10], a[13],
                            xbest, is_ae, use_rk4, *enables)
    obs = arrays[12]
    stats = calculate_stats(obs, est, best_logL, k)

    return {
        "species_key": species_key,
        "integrator": integrator,
        "config_key": config_key,
        "label": label,
        "full_name": full_name,
        "flags": {"washout_type": flags[0], "s_obs_estimated": flags[1],
                  "mortality_estimated": flags[2]},
        "temp_fit": bool(effective_temp_fit),
        "k": k,
        "active_parameters": [ALL_PARAM_NAMES[i] for i in active_idx],
        "xbest": xbest,
        "est": est,
        "obs": obs,
        "best_logL": best_logL,
        "best_nll": best_nll,
        "stats": stats,
        "elapsed": elapsed,
        "n_eval": int(n_eval_total),
        "scipy_success": bool(result.success),
        "scipy_message": str(result.message),
        "de_maxiter": int(de_maxiter),
        "de_popsize": int(de_popsize),
    }


# ============================================================
# OUTPUT
# ============================================================

# ============================================================
# MCMC ESTIMATION (emcee) + PHYSIOLOGICAL PRIOR
# ============================================================
# The thermal re-estimation is the part most exposed to "the parameters are
# free, so the fit means nothing" criticism. The defence (agreed design) is to
# keep the MODEL flexible but make the ESTIMATION PROTOCOL strict:
#   * sample the full posterior with emcee (not a single point), so the lack of
#     identifiability between confounded parameters shows up as posterior width
#     / correlation rather than a spuriously precise point estimate;
#   * constrain the thermal knobs through a PHYSIOLOGICAL PRIOR expressed on
#     ecologically meaningful quantities (the per-stage optimum and the warm-
#     side lethal temperature), declared BEFORE fitting against the published
#     Culex anchors (optimum 20-28 C, near-complete mortality ~35 C);
#   * report convergence (acceptance fraction, autocorrelation time) and
#     credible intervals.

try:
    import emcee
    EMCEE_AVAILABLE = True
except Exception:
    EMCEE_AVAILABLE = False


def _stage_vertex_summaries(x):
    """From a full parameter vector, return per-stage (T_opt, warm_lethal_T) for
    larva, pupa, adult using the vertex knobs at indices 15..23. Pure-python
    (used by the prior and by reporting); mirrors the njit kernel's math."""
    out = {}
    defs = (
        ("larva", AE_LARVA_TOPT, AE_LARVA_K, AE_LARVA_MMIN, M1_TW_MAX, 15),
        ("pupa",  AE_PUPA_TOPT,  AE_PUPA_K,  AE_PUPA_MMIN,  M1_TW_MAX, 18),
        ("adult", AE_ADULT_TOPT, AE_ADULT_K, AE_ADULT_MMIN, M3_TA_MAX, 21),
    )
    for name, ae_topt, ae_k, mmin, hi, base in defs:
        dTopt = float(x[base]); fwarm = float(x[base + 1])
        topt = ae_topt + dTopt
        kwarm = ae_k * fwarm
        if kwarm > 0:
            t_leth = topt + math.sqrt((1.0 / mmin) / kwarm)
        else:
            t_leth = float("inf")
        t_leth = min(t_leth, hi)
        out[name] = (topt, t_leth)
    return out


def log_prior_thermal(x, active_idx, bounds_full, temp_fit):
    """Log prior. Uniform (0) inside the box bounds for every active parameter,
    -inf outside. When thermal fitting is on, ADD a weak Gaussian penalty that
    pulls the per-stage optimum and the adult warm-lethal temperature toward the
    published Culex anchors. The penalty is deliberately weak (wide sigma): it
    keeps the sampler in physiological territory without overriding the data."""
    # Box bounds.
    for idx in active_idx:
        lo, hi = bounds_full[idx, 0], bounds_full[idx, 1]
        if not (lo <= x[idx] <= hi):
            return -np.inf
    if not temp_fit:
        return 0.0
    # Hard physiological admissibility on the optima of the ESTIMATED stages.
    summ = _stage_vertex_summaries(x)
    stage_for_base = {15: "larva", 18: "pupa", 21: "adult"}
    estimated_stage_names = {stage_for_base[b] for b in (15, 18, 21)
                             if b in cx_temp_active_idx()}
    for name in estimated_stage_names:
        topt, _ = summ[name]
        if not (CX_TOPT_HARD_LO <= topt <= CX_TOPT_HARD_HI):
            return -np.inf
    # Weak anchors (Culex): adult optimum centred in 20-28 band; adult warm-
    # lethal near 35 C. Wide sigmas => data-dominated.
    lp = 0.0
    ad_topt, ad_leth = summ["adult"]
    centre_opt = 0.5 * (CX_TOPT_LIT_LO + CX_TOPT_LIT_HI)   # 24 C
    sigma_opt = 0.5 * (CX_TOPT_LIT_HI - CX_TOPT_LIT_LO)    # 4 C (covers the band at 1 sigma)
    lp += -0.5 * ((ad_topt - centre_opt) / sigma_opt) ** 2
    sigma_leth = 3.0
    lp += -0.5 * ((ad_leth - CX_LETHAL_LIT) / sigma_leth) ** 2
    return lp


def run_one_mcmc(species_key, integrator, config_key, arrays,
                 temp_fit=False, n_walkers=None, n_steps=4000, n_burn=1000,
                 seed=123, progress=True):
    """Estimate one configuration by MCMC (emcee) and return a result dict in the
    same shape as run_one, plus posterior summaries. The point estimate reported
    is the posterior MEDIAN; credible intervals are the 16/84 percentiles."""
    if not EMCEE_AVAILABLE:
        raise ImportError("emcee is required for --mcmc (pip install emcee).")
    info = SPECIES_INFO[species_key]
    is_ae = (species_key == "Ae")
    use_rk4 = (integrator == "R")
    flags = CONFIG_FLAGS[config_key]
    bounds_full = info["bounds"]
    effective_temp_fit = bool(temp_fit) and (species_key == "Cx")
    global FIX_PHOTOPERIOD_SPECIES
    FIX_PHOTOPERIOD_SPECIES = species_key

    enable_diapause = 1
    enable_ccl = 1
    wash_code = WASH_CODE[flags[0]]
    enables = (enable_diapause, enable_ccl, wash_code)

    active_idx = get_active_indices(flags, effective_temp_fit)
    fixed_x = make_fixed_x(bounds_full, flags, effective_temp_fit)
    active_bounds = [(float(bounds_full[i, 0]), float(bounds_full[i, 1])) for i in active_idx]
    k = len(active_idx)

    label = build_label(config_key, integrator)
    if effective_temp_fit:
        label = label + ("+Ta" if TEMP_SCOPE == "adult" else "+T")
    label = label + "[MCMC]"
    full_name = build_full_name(config_key, integrator, info["diapause"])
    if effective_temp_fit:
        full_name = full_name + "; Culex thermal curves re-estimated (vertex form); MCMC (emcee)"

    a = arrays
    obs = a[12]

    def log_prob(theta):
        x = fixed_x.copy()
        for j, idx in enumerate(active_idx):
            x[idx] = theta[j]
        lp = log_prior_thermal(x, active_idx, bounds_full, effective_temp_fit)
        if not np.isfinite(lp):
            return -np.inf
        # b1 >= b2 diapause ordering (same constraint the likelihood enforces).
        if x[0] < x[1]:
            return -np.inf
        try:
            ll = float(simulate_loglik_direct(
                a[0], a[1], a[2], a[3], a[5], a[4],
                a[6], a[7], a[8], a[9], a[11], a[10], a[13], a[12],
                x, is_ae, use_rk4, *enables))
        except Exception:
            return -np.inf
        if not math.isfinite(ll):
            return -np.inf
        return lp + ll

    ndim = k
    if n_walkers is None:
        n_walkers = max(2 * ndim + 2, 24)
    rng = np.random.default_rng(seed)

    # Initialise walkers in a small ball around the bound midpoints (valid region).
    p0 = np.empty((n_walkers, ndim))
    mid = np.array([(lo + hi) / 2.0 for lo, hi in active_bounds])
    span = np.array([(hi - lo) for lo, hi in active_bounds])
    for w in range(n_walkers):
        ok = False
        while not ok:
            cand = mid + 0.05 * span * rng.standard_normal(ndim)
            cand = np.clip(cand,
                           [lo for lo, hi in active_bounds],
                           [hi for lo, hi in active_bounds])
            if np.isfinite(log_prob(cand)):
                p0[w] = cand
                ok = True

    print("\n" + "=" * 72)
    print(f"MCMC: {label}  ({info['pretty']})")
    print(f"  {full_name}")
    print(f"  walkers={n_walkers}  steps={n_steps}  burn={n_burn}  k={k}")
    print("=" * 72, flush=True)

    sampler = emcee.EnsembleSampler(n_walkers, ndim, log_prob)
    t0 = time.time()
    sampler.run_mcmc(p0, n_steps, progress=progress)
    elapsed = time.time() - t0

    # Convergence diagnostics.
    try:
        tau = sampler.get_autocorr_time(tol=0)
        tau_mean = float(np.nanmean(tau))
    except Exception:
        tau = None
        tau_mean = float("nan")
    acc = float(np.mean(sampler.acceptance_fraction))

    flat = sampler.get_chain(discard=n_burn, flat=True)
    logp = sampler.get_log_prob(discard=n_burn, flat=True)

    # Posterior summaries for active params.
    med = np.median(flat, axis=0)
    p16 = np.percentile(flat, 16, axis=0)
    p84 = np.percentile(flat, 84, axis=0)

    # Best (max-posterior) sample -> point vector for plotting / stats.
    ibest = int(np.argmax(logp))
    theta_best = flat[ibest]
    xbest = fixed_x.copy()
    for j, idx in enumerate(active_idx):
        xbest[idx] = theta_best[j]

    # Median vector (the reported estimate).
    xmed = fixed_x.copy()
    for j, idx in enumerate(active_idx):
        xmed[idx] = med[j]

    est = simulate_est_week(a[0], a[1], a[2], a[3], a[5], a[4],
                            a[6], a[7], a[8], a[9], a[11], a[10], a[13],
                            xmed, is_ae, use_rk4, *enables)
    # logL at the median (without prior) for fit statistics.
    best_logL = float(simulate_loglik_direct(
        a[0], a[1], a[2], a[3], a[5], a[4],
        a[6], a[7], a[8], a[9], a[11], a[10], a[13], a[12],
        xmed, is_ae, use_rk4, *enables))
    stats = calculate_stats(obs, est, best_logL, k)

    # Per-active credible intervals (raw scale).
    ci = {}
    for j, idx in enumerate(active_idx):
        ci[ALL_PARAM_NAMES[idx]] = (float(p16[j]), float(med[j]), float(p84[j]))

    print(f"\n  MCMC done in {elapsed:.1f}s")
    print(f"  acceptance fraction = {acc:.3f}  (target ~0.2-0.5)")
    print(f"  mean autocorr time  = {tau_mean:.1f} steps "
          f"(want n_steps >> ~50x this)")
    print(f"  posterior logL at median = {best_logL:.3f}")

    return {
        "species_key": species_key,
        "integrator": integrator,
        "config_key": config_key,
        "label": label,
        "full_name": full_name,
        "flags": {"washout_type": flags[0], "s_obs_estimated": flags[1],
                  "mortality_estimated": flags[2]},
        "temp_fit": bool(effective_temp_fit),
        "k": k,
        "active_parameters": [ALL_PARAM_NAMES[i] for i in active_idx],
        "xbest": xmed,                 # report the posterior median as the estimate
        "xmap": xbest,                 # also keep the max-a-posteriori vector
        "est": est,
        "obs": obs,
        "best_logL": best_logL,
        "best_nll": -best_logL,
        "stats": stats,
        "elapsed": elapsed,
        "n_eval": int(n_walkers * n_steps),
        "scipy_success": True,
        "scipy_message": "MCMC (emcee)",
        "de_maxiter": 0,
        "de_popsize": 0,
        "mcmc": {
            "n_walkers": n_walkers, "n_steps": n_steps, "n_burn": n_burn,
            "acceptance_fraction": acc, "autocorr_time_mean": tau_mean,
            "credible_intervals_16_50_84": ci,
        },
    }


# ============================================================
# Ae vs Cx STAGE MORTALITY CURVE PLOT
# ============================================================

def plot_ae_vs_cx_curves(res, out_dir):
    """Plot Aedes (original Jia) vs Culex (re-estimated) stage mortality curves
    for larva, pupa and adult, in one 3-panel figure. Uses the estimated vertex
    knobs in res['xbest']. Saved as PNG and PDF (publication-ready)."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        print(f"  (skipping curve plot: matplotlib unavailable: {e})")
        return None

    x = res["xbest"]
    defs = (
        ("Larva (m_L)", AE_LARVA_TOPT, AE_LARVA_MMIN, AE_LARVA_K,
         M1_TW_MIN, M1_TW_MAX, 15, "water"),
        ("Pupa (m_P)",  AE_PUPA_TOPT,  AE_PUPA_MMIN,  AE_PUPA_K,
         M2_TW_MIN, M1_TW_MAX, 18, "water"),
        ("Adult (m_A)", AE_ADULT_TOPT, AE_ADULT_MMIN, AE_ADULT_K,
         M3_TA_MIN, M3_TA_MAX, 21, "air"),
    )

    def curve(T, topt, mmin, kc, kw, lo, hi, cap=CX_MORT_CAP):
        T = np.asarray(T, float)
        k = np.where(T < topt, kc, kw)
        denom = (1.0 / mmin) - k * (T - topt) ** 2
        m = np.where(denom > 1e-12, 1.0 / denom, cap)
        m = np.clip(m, 0.0, cap)
        m = np.where((T > lo) & (T < hi), m, cap)
        return m

    Tg = np.linspace(0, 40, 500)
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
    for ax, (name, ae_topt, mmin, ae_k, lo, hi, base, ttype) in zip(axes, defs):
        # Aedes original (neutral knobs).
        m_ae = curve(Tg, ae_topt, mmin, ae_k, ae_k, lo, hi)
        # Culex estimated.
        dTopt = float(x[base]); fwarm = float(x[base + 1]); fcool = float(x[base + 2])
        cx_topt = ae_topt + dTopt
        m_cx = curve(Tg, cx_topt, mmin, ae_k * fcool, ae_k * fwarm, lo, hi)
        ax.plot(Tg, m_ae, color="#C0392B", lw=2.2, label="Ae (Jia 2016)")
        ax.plot(Tg, m_cx, color="#2471A3", lw=2.2, ls="--",
                label="Cx (re-estimated)")
        ax.axvspan(CX_TOPT_LIT_LO, CX_TOPT_LIT_HI, color="#2ecc71", alpha=0.10,
                   label="Cx optimum (lit. 20-28C)")
        ax.axvline(CX_LETHAL_LIT, color="gray", ls=":", lw=1.2,
                   label="Cx ~lethal (lit. ~35C)")
        ax.set_title(f"{name} [{ttype} temp]\nAe Topt={ae_topt:.1f} -> Cx Topt={cx_topt:.1f} C")
        ax.set_xlabel(f"{ttype} temperature (C)")
        ax.set_ylabel("daily mortality /day")
        ax.set_ylim(0, 0.5)
        ax.grid(alpha=0.25)
        if base == 15:
            ax.legend(fontsize=7, loc="upper left")
    fig.suptitle(f"Stage mortality: Aedes (Jia) vs Culex (re-estimated) -- {res['label']}",
                 fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    png = os.path.join(out_dir, f"{res['label'].replace('[','_').replace(']','')}_AE_vs_CX_curves.png")
    pdf = png.replace(".png", ".pdf")
    fig.savefig(png, dpi=130)
    fig.savefig(pdf)
    plt.close(fig)
    print(f"  Saved curve plot : {png}")
    return png


def res_flags_tuple(res):
    f = res["flags"]
    return (f["washout_type"], f["s_obs_estimated"], f["mortality_estimated"])


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
    active_set = set(int(i) for i in get_active_indices(res_flags_tuple(res),
                                                        res.get("temp_fit", False)))

    # Build a tidy parameter table (raw value, converted value where applicable,
    # and whether the parameter was estimated or held fixed for this config).
    converted_map = {2: ("cc", 10.0 ** xbest[2]),
                     3: ("a1", 10.0 ** xbest[3]),
                     4: ("a2", 10.0 ** xbest[4]),
                     5: ("pf_mm", 10.0 ** xbest[5]),
                     6: ("pf_scalar", 10.0 ** xbest[6]),
                     7: ("s_obs", 10.0 ** xbest[7]),
                     8: ("m_mult_egg", 10.0 ** xbest[8]),
                     9: ("m_mult_larva", 10.0 ** xbest[9]),
                     10: ("m_mult_pupa", 10.0 ** xbest[10]),
                     11: ("m_mult_adult", 10.0 ** xbest[11]),
                     12: ("m_mult_diapause", 10.0 ** xbest[12]),
                     13: ("pf_max", 10.0 ** xbest[13]),
                     14: ("hill", 10.0 ** xbest[14]),
                     15: ("larva_dTopt", xbest[15]), 16: ("larva_fwarm", xbest[16]), 17: ("larva_fcool", xbest[17]),
                     18: ("pupa_dTopt", xbest[18]),  19: ("pupa_fwarm", xbest[19]),  20: ("pupa_fcool", xbest[20]),
                     21: ("adult_dTopt", xbest[21]), 22: ("adult_fwarm", xbest[22]), 23: ("adult_fcool", xbest[23])}
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
        "optimizer": ("MCMC (emcee)" if res.get("mcmc") else "Differential Evolution"),
        "config_key": res["config_key"],
        "flags": res["flags"],
        "mcmc": res.get("mcmc"),
        "day_length": DAY_LENGTH,
        "delta_t": DELTA_T,
        "DE_MAXITER": res["de_maxiter"],
        "DE_POPSIZE": res["de_popsize"],
        "DE_SEEDS": list(DE_SEEDS),
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
            "pf_mm": float(10.0 ** xbest[5]),
            "pf_scalar": float(10.0 ** xbest[6]),
            "s_obs": float(10.0 ** xbest[7]),
            "dis": float(10.0 ** xbest[7]),
            "m_mult_egg": float(10.0 ** xbest[8]),
            "m_mult_larva": float(10.0 ** xbest[9]),
            "m_mult_pupa": float(10.0 ** xbest[10]),
            "m_mult_adult": float(10.0 ** xbest[11]),
            "m_mult_diapause": float(10.0 ** xbest[12]),
            "pf_max": float(10.0 ** xbest[13]),
            "hill": float(10.0 ** xbest[14]),
            "larva_dTopt": float(xbest[15]), "larva_fwarm": float(xbest[16]), "larva_fcool": float(xbest[17]),
            "pupa_dTopt": float(xbest[18]),  "pupa_fwarm": float(xbest[19]),  "pupa_fcool": float(xbest[20]),
            "adult_dTopt": float(xbest[21]), "adult_fwarm": float(xbest[22]), "adult_fcool": float(xbest[23]),
        },
        "temp_fit": bool(res.get("temp_fit", False)),
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
    print(f"    pf (mm)   = {10.0 ** xbest[5]:.6g}   (washout threshold / Hill half-saturation, mm of driver)")
    print(f"    pf_scalar = {10.0 ** xbest[6]:.6g}   (linear washout only)")
    print(f"    pf_max    = {10.0 ** xbest[13]:.6g}   (Hill washout max rate/day)")
    print(f"    hill (h)  = {10.0 ** xbest[14]:.6g}   (Hill exponent)")
    print(f"    s_obs     = {10.0 ** xbest[7]:.6g}   (= dis)")
    print(f"    m_mult    = egg {10.0**xbest[8]:.4g}, larva {10.0**xbest[9]:.4g}, "
          f"pupa {10.0**xbest[10]:.4g}, adult {10.0**xbest[11]:.4g}, "
          f"diap {10.0**xbest[12]:.4g}")

    if res.get("temp_fit", False):
        print("\n  Culex thermal curves (vertex re-parameterisation):")
        stage_defs = (
            ("larva", AE_LARVA_TOPT, AE_LARVA_K, M1_TW_MIN, M1_TW_MAX, 15),
            ("pupa",  AE_PUPA_TOPT,  AE_PUPA_K,  M2_TW_MIN, M1_TW_MAX, 18),
            ("adult", AE_ADULT_TOPT, AE_ADULT_K, M3_TA_MIN, M3_TA_MAX, 21),
        )
        all_ok = True
        for sname, ae_topt, ae_k, lo, hi, base in stage_defs:
            dTopt = float(xbest[base]); fwarm = float(xbest[base+1]); fcool = float(xbest[base+2])
            topt = ae_topt + dTopt
            kwarm = ae_k * fwarm
            mmin = (AE_LARVA_MMIN if sname=="larva" else
                    AE_PUPA_MMIN if sname=="pupa" else AE_ADULT_MMIN)
            # Warm-side lethal temperature: where vertex denom hits zero ->
            # T_lethal = T_opt + sqrt( (1/m_min) / k_warm ).
            import math as _m
            t_lethal = topt + _m.sqrt((1.0/mmin)/kwarm) if kwarm > 0 else float("inf")
            t_lethal = min(t_lethal, hi)   # window cap also forces lethal
            print(f"    {sname:5s}: T_opt {ae_topt:.1f}->{topt:.1f} C  "
                  f"k_warm x{fwarm:.2f}  k_cool x{fcool:.2f}  "
                  f"warm-lethal ~{t_lethal:.1f} C")
        # Literature validation against Cx anchors (declared a priori).
        ad_topt = AE_ADULT_TOPT + float(xbest[21])
        opt_ok = (CX_TOPT_LIT_LO <= ad_topt <= CX_TOPT_LIT_HI)
        # adult warm-lethal
        ad_kw = AE_ADULT_K * float(xbest[22])
        import math as _m
        ad_lethal = AE_ADULT_TOPT + float(xbest[21]) + _m.sqrt((1.0/AE_ADULT_MMIN)/ad_kw) if ad_kw>0 else float("inf")
        ad_lethal = min(ad_lethal, M3_TA_MAX)
        leth_ok = (CX_LETHAL_LIT - 3.0 <= ad_lethal <= CX_LETHAL_LIT + 3.0)
        print("    --- literature check (declared before fitting) ---")
        print(f"    adult T_opt in {CX_TOPT_LIT_LO:.0f}-{CX_TOPT_LIT_HI:.0f} C : "
              f"{'OK' if opt_ok else 'OUTSIDE'}  (got {ad_topt:.1f} C)")
        print(f"    adult warm-lethal near ~{CX_LETHAL_LIT:.0f} C : "
              f"{'OK' if leth_ok else 'OUTSIDE'}  (got {ad_lethal:.1f} C)")
        if not (opt_ok and leth_ok):
            print("    NOTE: a value outside the band is not necessarily wrong -- it may")
            print("          reflect Japanese-population thermal adaptation rather than a")
            print("          model error. Report it; do not silently force it into range.")

    print(f"\n  Saved JSON       : {json_path}")
    print(f"  Saved weekly CSV : {csv_path}")
    print(f"  Saved param CSV  : {param_csv_path}")

    # Aedes-vs-Culex stage mortality curve plot (only meaningful when the Cx
    # thermal knobs were estimated). Controlled by a module flag set from CLI.
    if res.get("temp_fit", False) and not globals().get("_NO_CURVE_PLOT", False):
        try:
            plot_ae_vs_cx_curves(res, out_dir)
        except Exception as e:
            print(f"  (curve plot failed: {e})")

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
            "washout_type": o["flags"]["washout_type"],
            "s_obs_estimated": o["flags"]["s_obs_estimated"],
            "mortality_estimated": o["flags"]["mortality_estimated"],
            "temp_fit": o.get("temp_fit", False),
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
    any_temp = any(o.get("temp_fit", False) for o in all_outs)
    suffix = "_T" if any_temp else ""
    out_csv = os.path.join(script_dir, f"{info['code']}_M{integrator}{suffix}_summary.csv")
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
    parser.add_argument("--site", choices=list(SITE_CONFIG.keys()), default="tokyo",
                        help="data site: tokyo (2003-2013) or saitama (2014-2019, "
                             "Kumagaya climate). Default tokyo. Selects climate/obs "
                             "paths and year ranges; the model itself is unchanged.")
    parser.add_argument("--fix-photoperiod", action="store_true",
                        help="Hold the photoperiod parameters (b1, b2, a1, a2) at the "
                             "Tokyo-estimated values instead of re-estimating them. Use "
                             "when transferring the Tokyo diapause timing to Saitama so "
                             "that only landscape-dependent (and, with --temp-fit, thermal) "
                             "parameters are re-estimated.")
    parser.add_argument("--temp-fit", action="store_true",
                        help="(Cx only) additionally estimate the Culex vertex thermal "
                             "knobs. Equivalent to appending '+T' to the config key. The "
                             "number of knobs depends on --temp-stages. Ignored for --species Ae.")
    parser.add_argument("--temp-stages", choices=["all", "adult"], default="adult",
                        help="Which life stages get their thermal curve re-estimated when "
                             "--temp-fit is on. 'adult' (default) re-estimates only the adult "
                             "curve (3 knobs) and holds larva/pupa at the Aedes curves; this is "
                             "recommended because the adult thermal response is well identified "
                             "by the adult observations while the water-stage knobs are not. "
                             "'all' re-estimates larva+pupa+adult (9 knobs).")
    parser.add_argument("--mcmc", action="store_true",
                        help="Estimate by MCMC (emcee) instead of Differential Evolution. "
                             "Returns the full posterior, convergence diagnostics, and 16/50/84 "
                             "credible intervals. Recommended for the Culex thermal re-estimation.")
    parser.add_argument("--mcmc-steps", type=int, default=4000, help="MCMC steps per walker.")
    parser.add_argument("--mcmc-burn", type=int, default=1000, help="MCMC burn-in steps discarded.")
    parser.add_argument("--mcmc-walkers", type=int, default=0,
                        help="MCMC walkers (0 = auto: max(2k+2, 24)).")
    parser.add_argument("--no-curve-plot", action="store_true",
                        help="Skip the Aedes-vs-Culex stage mortality curve plot.")
    args = parser.parse_args()

    # Activate the site BEFORE any njit kernel compiles (WEEK_LENGTH is baked in).
    apply_site_config(args.site)

    global TEMP_SCOPE
    TEMP_SCOPE = args.temp_stages
    if TEMP_SCOPE == "adult":
        print("Thermal re-estimation scope: ADULT ONLY (larva & pupa held at the "
              "Aedes curves). The water-stage thermal knobs are not identified by "
              "the adult observations, so only the adult curve is re-estimated.",
              flush=True)
    else:
        print("Thermal re-estimation scope: ALL stages (larva + pupa + adult, 9 knobs). "
              "Note: the water-stage knobs are typically poorly identified.", flush=True)

    global FIX_PHOTOPERIOD
    FIX_PHOTOPERIOD = bool(args.fix_photoperiod)
    if FIX_PHOTOPERIOD:
        if args.site == "tokyo":
            print("WARNING: --fix-photoperiod has no effect at the Tokyo calibration "
                  "site; photoperiod is being held at the Tokyo values it was "
                  "estimated from.", flush=True)
        print("Photoperiod parameters (b1, b2, a1, a2) are FIXED at the Tokyo values; "
              "only landscape-dependent (and thermal, if --temp-fit) parameters will "
              "be estimated.", flush=True)

    species_key = args.species or prompt_choice("Species", ["Ae", "Cx"])
    integrator = args.integrator or prompt_choice("Integrator", ["E", "R"])
    if args.config:
        config_choice = args.config
    else:
        print("Single-model configs:", ", ".join(CONFIG_ORDER))
        print("  (MA = full all-components model)")
        print("  (append '+T' to estimate the Culex thermal response, e.g. MM+T)")
        print("Batch mode: 'all'  -> estimate every config above in turn + summary CSV")
        print("            'all+T' -> same batch, Culex thermal response estimated")
        print("            (this is a run mode, not a model; MA is included once)")
        config_choice = input("Config or 'all': ").strip()

    # A trailing '+T' on the config selects Culex thermal re-estimation; it is
    # equivalent to --temp-fit. Strip it off to recover the base config key.
    temp_fit = bool(args.temp_fit)
    if config_choice.endswith("+T"):
        temp_fit = True
        config_choice = config_choice[:-2]
    if temp_fit and species_key != "Cx":
        print("Note: --temp-fit / +T has no effect for species Ae; ignoring it.")
        temp_fit = False

    info = SPECIES_INFO[species_key]
    print(f"\nLoading input data for {info['pretty']} at {POINT} ...")
    arrays = load_input_arrays(info["full"], info["obs_candidates"], species_code=info["code"])
    obs = arrays[12]
    if np.sum(obs >= 0.0) == 0:
        raise RuntimeError("No valid observation data found.")
    print(f"Valid observation weeks: {int(np.sum(obs >= 0.0))}")

    # Curve-plot toggle (read by write_outputs).
    globals()["_NO_CURVE_PLOT"] = bool(args.no_curve_plot)
    use_mcmc = bool(args.mcmc)
    mcmc_walkers = None if args.mcmc_walkers == 0 else int(args.mcmc_walkers)

    def estimate_one(ck):
        if use_mcmc:
            return run_one_mcmc(species_key, integrator, ck, arrays,
                                temp_fit=temp_fit, n_walkers=mcmc_walkers,
                                n_steps=args.mcmc_steps, n_burn=args.mcmc_burn)
        return run_one(species_key, integrator, ck, arrays,
                       args.de_maxiter, args.de_popsize, temp_fit)

    if config_choice.lower() == "all":
        all_outs = []
        for ck in CONFIG_ORDER:
            res = estimate_one(ck)
            _, _, out = write_outputs(res)
            all_outs.append(out)
        write_summary(all_outs, species_key, integrator)
    else:
        if config_choice not in CONFIG_FLAGS:
            raise SystemExit(f"Unknown config '{config_choice}'. "
                             f"Choose from: {', '.join(CONFIG_ORDER)} or 'all' "
                             f"(optionally with a '+T' suffix for Cx thermal fitting).")
        res = estimate_one(config_choice)
        write_outputs(res)

    print("\nDone.")


if __name__ == "__main__":
    main()