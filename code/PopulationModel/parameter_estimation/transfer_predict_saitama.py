#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
transfer_predict_saitama.py

東京で較正したパラメータ（転移）と、埼玉で再較正したパラメータ（再推定）の
両方を読み込み、埼玉の気候データに適用して予測系列を計算、実測と並べた図を描くスクリプト。

このスクリプトは推定（DE最適化）を一切行わない。やることは:
  1) 推定済みの JSON から 15 個の生パラメータ x を読み込む（東京用・埼玉用の両方）
  2) 埼玉の気候データを使って main.py と同一のモデルを前進積分する
  3) 週ごとの成虫予測 est_week を作る
  4) 実測と R^2_log / RMSE_log / r を計算し、比較図にする
  ※ 視認性向上のため、Mモデルは完全に除外しています。
"""

import os
import sys
import json
import math
import argparse

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ============================================================
# 設定（必要に応じて書き換える）
# ============================================================

# --- 埼玉データの場所と期間 -------------------------------------------------
DATA_BASE       = "./../.."         # input_data / measurement_data の親
POINT           = "kumagaya"
OBSPOINT        = "saitama"
START_YEAR      = 2014              # 埼玉の評価開始年
YEAR_LENGTH     = 6                 # 2014..2019 の 6 年
WARMUP_START    = 2002              # ウォームアップ開始年
WARMUP_END      = 2013              # ウォームアップ終了年

# --- 推定パラメータ JSON の場所 -----------------------------------------
PARAM_DIR       = "."               # *_params.json がある場所
INTEGRATOR      = "R"               # 図に使う積分（R=RK4 / E=Euler）

# JSON ファイル名（種 × モデル）。None のままなら PARAM_DIR から自動探索する。
# "_S" がつくものは埼玉で再較正されたパラメータ用
PARAM_FILES = {
    ("Ae", "B_PCMP"): None, ("Ae", "MM"): None,
    ("Ae", "MM_S"): None,
    ("Cx", "B_PCMP"): None, ("Cx", "MM"): None,
    ("Cx", "MM_S"): None,
}

# ============================================================
# 以降は main.py から移植した定数・関数
# ============================================================

DAY_LENGTH = 50
DELTA_T = 1.0 / DAY_LENGTH

# 物理定数
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
M4_BASE_CX = 0.005

WASH_CODE = {"off": 0, "threshold": 1, "linear": 2, "hill": 3}

ALL_PARAM_NAMES = [
    "p_b1", "p_b2", "p_ccl", "p_a1l", "p_a2l", "p_pf",
    "p_pf_scalar_l", "p_dis_l",
    "p_m_mult_egg_l", "p_m_mult_larva_l", "p_m_mult_pupa_l",
    "p_m_mult_adult_l", "p_m_mult_diapause_l",
    "p_pf_max_l", "p_hill_l",
]

SPECIES_INFO = {
    "Ae": {"code": "Ae", "full": "Aedes_albopictus", "pretty": "Aedes albopictus",
           "is_ae": True,
           "obs_candidates": ["Aealbopictus_female", "Aealbo_female",
                              "Aedes_albopictus_female", "Aedes_albopictus", "Aealbo"]},
    "Cx": {"code": "Cx", "full": "Culex_pipiens", "pretty": "Culex pipiens",
           "is_ae": False,
           "obs_candidates": ["Cxpipiens_female", "Culex_pipiens_female",
                              "Cx_pipiens_female", "Cxpipiens", "Culex_pipiens"]},
}


# ---- 生理関数 ----------------------------------------------------------------

def safe_exp_arg(x):
    if x > 700.0: return 700.0
    if x < -700.0: return -700.0
    return x

def clamp01_08(x):
    if x < 0.0: return 0.0
    if x > 0.8: return 0.8
    return x

def cal_diapause(dwk, a1, b1, a2, b2, enable):
    if enable == 0: return 0.0, 0.0
    z0 = 1.0 / (1.0 + math.exp(safe_exp_arg(-a1 * (b1 - dwk))))
    z1 = 1.0 / (1.0 + math.exp(safe_exp_arg(a2 * (b2 - dwk))))
    return z0, z1

def cal_carrying_capacity(pn, cc, enable):
    ccl = cc * pn if enable else cc
    if ccl <= 0.0: ccl = 1.0e-6
    return ccl

def cal_development_rate(ta, tw):
    d0 = max(0.0, 0.507 * math.exp(-1.0 * ((tw - 30.85) / 12.82))) if tw > D0_TW_MIN else 1.0 / 60.0
    d1 = max(0.0, 0.1727 * math.exp(-1.0 * ((tw - 28.40) / 10.20))) if tw > D1_TW_MIN else 1.0 / 60.0
    d2 = max(0.0, 0.602 * math.exp(-1.0 * ((tw - 34.29) / 15.07))) if tw > D2_TW_MIN else 1.0 / 60.0
    d3 = max(0.0, -15.837 + 1.2897 * ta - 0.0163 * (ta * ta)) if ta > D3_TA_MIN else 0.0
    return d0, d1, d2, d3

def cal_mortality_rate(ta, tw):
    m0 = 0.05 if M0_TW_MIN < tw < 31.12 else 1.0
    if M1_TW_MIN < tw < M1_TW_MAX:
        denom1 = -0.1305 * (tw * tw) + 3.868 * tw + 30.83
        denom2 = -0.1502 * (tw * tw) + 5.057 * tw + 3.517
        m1 = max(0.0, 1.0 / denom1 if abs(denom1) > 1.0e-12 else 1.0)
        m2 = max(0.0, 1.0 / denom2 if abs(denom2) > 1.0e-12 else 1.0)
    else:
        m1, m2 = 1.0, 1.0
    if M3_TA_MIN < ta < M3_TA_MAX:
        denom3 = -0.1921 * (ta * ta) + 8.147 * ta - 22.98
        m3 = max(0.0, 1.0 / denom3 if abs(denom3) > 1.0e-12 else 1.0)
    else:
        m3 = 1.0
    return m0, m1, m2, m3

def cal_washout_mortality(drive, pf, pf_scalar, pf_max, hill, wash_code):
    if wash_code == 0: return 0.0
    if wash_code == 1: return 0.0 if drive <= pf else 1.0
    if wash_code == 2:
        if drive <= pf: return 0.0
        return min(0.8, max(0.0, pf_scalar * (drive - pf)))
    # hill
    if drive <= 0.0: return 0.0
    denom = (drive ** hill) + (pf ** hill)
    if denom <= 0.0: return 0.0
    return min(0.95, max(0.0, pf_max * (drive ** hill) / denom))

def unpack_params(x):
    b1, b2 = x[0], x[1]
    cc = 10.0 ** x[2]
    a1, a2 = 10.0 ** x[3], 10.0 ** x[4]
    pf = 10.0 ** x[5]
    pf_scalar = 10.0 ** x[6]
    dis = 10.0 ** x[7]
    m_mult_egg = 10.0 ** x[8]
    m_mult_larva = 10.0 ** x[9]
    m_mult_pupa = 10.0 ** x[10]
    m_mult_adult = 10.0 ** x[11]
    m_mult_diapause = 10.0 ** x[12]
    pf_max = 10.0 ** x[13]
    hill = 10.0 ** x[14]
    return (b1, b2, a1, a2, cc, pf, pf_scalar, dis,
            m_mult_egg, m_mult_larva, m_mult_pupa, m_mult_adult, m_mult_diapause,
            pf_max, hill)

def dynamics(y, is_ae, z0, z1, ccl, d0, d1, d2, d3, m0, m1, m2, m3, m4, mpf):
    x0, x1, x2, x3, x4 = y[0], y[1], y[2], y[3], y[4]
    dy = np.zeros(5, dtype=np.float64)
    if is_ae:
        dy[0] = (1.0 - z0) * d3 * x3 - (m0 + d0 + mpf) * x0 + z1 * x4
        dy[1] = d0 * x0 - ((m1 + d1 + mpf) + (x1 / ccl)) * x1
        dy[2] = d1 * x1 - (m2 + d2 + mpf) * x2
        dy[3] = d2 * x2 - m3 * x3
        dy[4] = z0 * d3 * x3 - (z1 + m4) * x4
    else:
        dy[0] = d3 * x3 - (m0 + d0 + mpf) * x0
        dy[1] = d0 * x0 - ((m1 + d1 + mpf) + (x1 / ccl)) * x1
        dy[2] = d1 * x1 - (m2 + d2 + mpf) * x2
        dy[3] = d2 * x2 - (m3 + z0) * x3 + z1 * x4
        dy[4] = z0 * x3 - (z1 + m4) * x4
    return dy

def step_once(stage_n, is_ae, use_rk4, z0, z1, ccl,
              d0, d1, d2, d3, m0, m1, m2, m3, m4, mpf):
    if use_rk4:
        k1 = dynamics(stage_n, is_ae, z0, z1, ccl, d0, d1, d2, d3, m0, m1, m2, m3, m4, mpf)
        y2 = stage_n + 0.5 * DELTA_T * k1
        k2 = dynamics(y2, is_ae, z0, z1, ccl, d0, d1, d2, d3, m0, m1, m2, m3, m4, mpf)
        y3 = stage_n + 0.5 * DELTA_T * k2
        k3 = dynamics(y3, is_ae, z0, z1, ccl, d0, d1, d2, d3, m0, m1, m2, m3, m4, mpf)
        y4 = stage_n + DELTA_T * k3
        k4 = dynamics(y4, is_ae, z0, z1, ccl, d0, d1, d2, d3, m0, m1, m2, m3, m4, mpf)
        stage_n += (DELTA_T / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
    else:
        k1 = dynamics(stage_n, is_ae, z0, z1, ccl, d0, d1, d2, d3, m0, m1, m2, m3, m4, mpf)
        stage_n += DELTA_T * k1
    for i in range(5):
        if stage_n[i] < 0.0: stage_n[i] = 0.0


# ============================================================
# データローダー
# ============================================================

def _read_climate_6(fname):
    df = pd.read_csv(fname, header=None)
    df = df.iloc[:, :6]
    df.columns = ["TA", "TW", "DWK", "SMOI", "PREC", "ROFF"]
    return df

def load_input_arrays(species_code, species_full, obs_candidates, data_base, point, obs_point,
                      start_year, year_length, warmup_start, warmup_end):
    ta_w, tw_w, dwk_w, smoi_w, prec_w, roff_w = [], [], [], [], [], []
    for yr in range(warmup_start, warmup_end + 1):
        cands = [f"{data_base}/input_data/{point}/{yr}/{point}{yr}_climdata2.csv",
                 f"{data_base}/input_data/{point}/{yr}/{point}{yr}_climdata.csv"]
        fname = next((c for c in cands if os.path.exists(c)), None)
        if fname is None: continue
        df = _read_climate_6(fname)
        ta_w += df["TA"].tolist(); tw_w += df["TW"].tolist()
        dwk_w += df["DWK"].tolist(); smoi_w += df["SMOI"].tolist()
        prec_w += df["PREC"].tolist(); roff_w += df["ROFF"].tolist()

    ta, tw, dwk, smoi, prec, roff = [], [], [], [], [], []
    data_year = [0] * year_length
    for yr in range(start_year, start_year + year_length):
        cands = [f"{data_base}/input_data/{point}/{yr}/{point}{yr}_climdata.csv",
                 f"{data_base}/input_data/{point}/{yr}/{point}{yr}_climdata2.csv"]
        fname = next((c for c in cands if os.path.exists(c)), None)
        if fname is None: continue
        df = _read_climate_6(fname)
        data_year[yr - start_year] = len(df)
        ta += df["TA"].tolist(); tw += df["TW"].tolist()
        dwk += df["DWK"].tolist(); smoi += df["SMOI"].tolist()
        prec += df["PREC"].tolist(); roff += df["ROFF"].tolist()

    week_length = 52 * year_length
    obs_week = np.full(week_length, -999.0, dtype=np.float64)
    
    week_plus = 0
    obs_found_any = False
    base_dir = f"{data_base}/measurement_data/{obs_point}/{species_full}"

    for yr in range(start_year, start_year + year_length):
        if yr != start_year: week_plus += 52
            
        cands = [f"{base_dir}/{species_code}_{yr}_{obs_point}_obs.txt"]
        for s1 in obs_candidates:
            cands.append(f"{base_dir}/{obs_point}_{s1}_{yr}.tsv")
            
        fname = next((c for c in cands if os.path.exists(c)), None)
        
        if fname:
            if not obs_found_any:
                print(f"Obs file found: {fname}")
                obs_found_any = True
                
            df = pd.read_csv(fname, sep="\t", header=None,
                             names=["week", "pop"], engine="python")
            for _, row in df.iterrows():
                wc = int(row["week"]) + week_plus - 1
                if 0 <= wc < week_length:
                    obs_week[wc] = float(row["pop"])
        else:
            print(f"Warning: observation file missing for {yr} in {base_dir}")

    npf = lambda v: np.asarray(v, dtype=np.float64)
    return {
        "wup": (npf(ta_w), npf(tw_w), npf(dwk_w), npf(smoi_w), npf(roff_w), npf(prec_w)),
        "main": (npf(ta), npf(tw), npf(dwk), npf(smoi), npf(roff), npf(prec)),
        "obs": obs_week,
        "data_year": np.asarray(data_year, dtype=np.int64),
        "week_length": week_length,
    }


# ============================================================
# 予測
# ============================================================

def simulate_est_week(arrays, x, is_ae, use_rk4, enable_diapause, enable_ccl, wash_code):
    (b1, b2, a1, a2, cc, pf, pf_scalar, dis,
     m_egg, m_lar, m_pup, m_adt, m_dia, pf_max, hill) = unpack_params(x)
    week_length = arrays["week_length"]

    if enable_diapause == 1 and b1 < b2:
        return np.zeros(week_length, dtype=np.float64)

    stage_n = np.zeros(5, dtype=np.float64)
    stage_n[3] = 100.0
    stage_n[4] = 100.0

    def mortalities(ta_v, tw_v):
        m0_raw, m1, m2, m3 = cal_mortality_rate(ta_v, tw_v)
        m0 = clamp01_08(m0_raw * m_egg)
        m1 = clamp01_08(m1 * m_lar)
        m2 = clamp01_08(m2 * m_pup)
        m3 = clamp01_08(m3 * m_adt)
        m4 = clamp01_08(m0_raw * m_dia) if is_ae else clamp01_08(M4_BASE_CX * m_dia)
        return m0, m1, m2, m3, m4

    # ---- warmup ----
    ta_w, tw_w, dwk_w, smoi_w, roff_w, prec_w = arrays["wup"]
    for day in range(ta_w.shape[0]):
        z0, z1 = cal_diapause(dwk_w[day], a1, b1, a2, b2, enable_diapause)
        ccl = cal_carrying_capacity(smoi_w[day] / 100.0, cc, enable_ccl)
        d0, d1, d2, d3 = cal_development_rate(ta_w[day], tw_w[day])
        m0, m1, m2, m3, m4 = mortalities(ta_w[day], tw_w[day])
        mpf = cal_washout_mortality(roff_w[day], pf, pf_scalar, pf_max, hill, wash_code)
        for _ in range(DAY_LENGTH):
            step_once(stage_n, is_ae, use_rk4, z0, z1, ccl,
                      d0, d1, d2, d3, m0, m1, m2, m3, m4, mpf)

    # ---- main ----
    ta, tw, dwk, smoi, roff, prec = arrays["main"]
    data_year = arrays["data_year"]
    est_week = np.zeros(week_length, dtype=np.float64)
    stage_sum = np.zeros(5, dtype=np.float64)
    y, day_count_for_weekave, week_count, year_day_count = 0, 1, 0, 0

    for day in range(ta.shape[0]):
        z0, z1 = cal_diapause(dwk[day], a1, b1, a2, b2, enable_diapause)
        ccl = cal_carrying_capacity(smoi[day] / 100.0, cc, enable_ccl)
        d0, d1, d2, d3 = cal_development_rate(ta[day], tw[day])
        m0, m1, m2, m3, m4 = mortalities(ta[day], tw[day])
        mpf = cal_washout_mortality(roff[day], pf, pf_scalar, pf_max, hill, wash_code)
        for _ in range(DAY_LENGTH):
            step_once(stage_n, is_ae, use_rk4, z0, z1, ccl, d0, d1, d2, d3, m0, m1, m2, m3, m4, mpf)
            stage_sum += stage_n
        
        if day_count_for_weekave % 7 == 0:
            if week_count < week_length:
                est_week[week_count] = dis * (stage_sum[3] / (DAY_LENGTH * 7.0))
            week_count += 1
            stage_sum[:] = 0.0
            
        day_count_for_weekave += 1
        year_day_count += 1
        
        if y < data_year.shape[0] and year_day_count >= data_year[y]:
            y += 1
            year_day_count = 0
            day_count_for_weekave = 1
            stage_sum[:] = 0.0
            if y >= data_year.shape[0]: break
        if week_count >= week_length: break
            
    return est_week


# ============================================================
# パラメータ JSON 読み込み
# ============================================================

def x_from_json(path):
    with open(path, "r", encoding="utf-8") as f:
        j = json.load(f)
    params = j["parameters"]
    x = np.array([float(params[name]) for name in ALL_PARAM_NAMES], dtype=np.float64)
    flags = j.get("flags", {})
    wt = flags.get("washout_type", None)
    if wt is None:
        lab = j.get("label", "")
        # B_PCMP and MM both use the original threshold washout; only the
        # washout-augmented configs (MA / MW / MA-*) use the Hill form.
        if "B_PCMP" in lab or lab.startswith(("MRM", "MEM")) or "MM" in lab:
            wt = "threshold"
        else:
            wt = "hill"
    wash_code = WASH_CODE.get(wt, 1)
    integrator = j.get("integration", "RK4")
    use_rk4 = (integrator.upper().startswith("RK") or integrator == "R")
    return x, wash_code, use_rk4, j.get("label", os.path.basename(path))

def autodiscover(param_dir, species_code, model_key, integrator, target_region="tokyo"):
    cand_labels = []
    if model_key == "B_PCMP":
        cand_labels = [f"B_PCMP_{integrator}", "B_PCMP_R", "B_PCMP_E", "B_PCMP"]
    elif model_key == "MM":
        cand_labels = [f"M{integrator}M", "MRM", "MEM"]
    elif model_key == "MA":
        cand_labels = [f"M{integrator}A", "MRA", "MEA", "MRA-washout", "MEA-washout"]

    # パターン1: 厳密な地域名付き
    for lab in cand_labels:
        tok = f"{species_code}_{lab}_{target_region}_params.json"
        p = os.path.join(param_dir, tok)
        if os.path.exists(p): return p

    # パターン2: 曖昧検索（ファイル名に地域名とラベル名が含まれているか）
    for fn in os.listdir(param_dir):
        if not fn.endswith(".json"): continue
        if fn.startswith(species_code + "_") and target_region in fn and any(lab in fn for lab in cand_labels):
            return os.path.join(param_dir, fn)

    return None


# ============================================================
# メトリクスと図
# ============================================================

def metrics(obs, est):
    obs, est = np.asarray(obs, dtype=float), np.asarray(est, dtype=float)
    valid = (obs >= 0.0) & np.isfinite(obs) & np.isfinite(est)
    n = int(np.sum(valid))
    if n == 0: return {"n": 0, "R2_Log10": float("nan"), "RMSE_Log10": float("nan"), "r": float("nan")}
    y = obs[valid]; yh = est[valid]
    ly, lyh = np.log10(y + 1.0), np.log10(yh + 1.0)
    ss_res = float(np.sum((ly - lyh) ** 2))
    ss_tot = float(np.sum((ly - np.mean(ly)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    rmse = float(np.sqrt(np.mean((ly - lyh) ** 2)))
    r = float(np.corrcoef(ly, lyh)[0, 1]) if np.std(ly) > 0 and np.std(lyh) > 0 else float("nan")
    return {"n": n, "R2_Log10": r2, "RMSE_Log10": rmse, "r": r}

# Figure shows: B_PCMP (legacy) transferred, MM (best model) transferred, and
# MM re-estimated at Saitama. MA is no longer the focus.
MODEL_STYLE = {
    "B_PCMP": dict(color="#1f4e79", label="B_PCMP (Tokyo transfer)", lw=1.6, ls="-"),
    "MM":     dict(color="#c0392b", label="MM (Tokyo transfer)", lw=1.8, ls="-"),
    "MM_S":   dict(color="#c0392b", label="MM (Saitama re-estimation)", lw=1.8, ls="--", alpha=0.85),
}
MODEL_ORDER = ["B_PCMP", "MM", "MM_S"]


def build_week_year_axis(data_year, start_year, week_length):
    xs = np.full(week_length, np.nan)
    wk = 0
    for yi, ndays in enumerate(data_year):
        for w in range(52):
            if wk >= week_length: break
            xs[wk] = (start_year + yi) + (w / 52.0)
            wk += 1
    return xs

def plot_timeseries(results, start_year, out_path):
    # 2 rows (species) x 3 cols (Plan B):
    #   col 0: B_PCMP transferred from Tokyo (legacy model, no re-fit)
    #   col 1: MM transferred from Tokyo (revised best model, no re-fit)
    #   col 2: MM re-estimated at Saitama (landscape parameters re-estimated)
    species = ["Ae", "Cx"]
    # (key in results["est"], colour, line style, column title)
    panels = [("B_PCMP", "#1f4e79", "-",  "B_PCMP\n(Tokyo transfer)"),
              ("MM",     "#c0392b", "-",  "MM\n(Tokyo transfer)"),
              ("MM_S",   "#c0392b", "--", "MM\n(Saitama re-estimation)")]
    fig, axes = plt.subplots(2, 3, figsize=(13.5, 7.2), sharex=True, sharey="row")
    for ri, sk in enumerate(species):
        info = SPECIES_INFO[sk]
        R = results[sk]
        xs = R["xaxis"]; obs = R["obs"]; valid = obs >= 0.0
        m = R["metrics"]
        for ci, (mk, color, ls, title) in enumerate(panels):
            ax = axes[ri, ci]
            ax.scatter(xs[valid], np.log10(obs[valid] + 1.0), s=7, color="0.5",
                       alpha=0.45, zorder=2, linewidths=0)
            if mk in R["est"]:
                ax.plot(xs, np.log10(R["est"][mk] + 1.0), color=color, ls=ls, lw=1.6, zorder=3)
            if mk in m:
                txt = f"R\u00b2={m[mk]['R2_Log10']:.2f}\nRMSE={m[mk]['RMSE_Log10']:.2f}"
                ax.text(0.985, 0.96, txt, transform=ax.transAxes, ha="right", va="top",
                        fontsize=7.6, family="monospace",
                        bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="0.8", lw=0.7))
            ax.grid(True, ls=":", alpha=0.4); ax.set_axisbelow(True)
            if ri == 0: ax.set_title(title, fontsize=10.5)
            if ci == 0:
                g, s = info['pretty'].split()[0], info['pretty'].split()[1]
                ax.set_ylabel(f"$\\it{{{g}}}$ $\\it{{{s}}}$\nlog10 (adult + 1)", fontsize=9.5)
            if ri == 1: ax.set_xlabel("Year", fontsize=9.5)
    from matplotlib.lines import Line2D
    handles = [Line2D([0],[0],color="#1f4e79",lw=2,label="B_PCMP (legacy)"),
               Line2D([0],[0],color="#c0392b",lw=2,label="MM (revised, best model)"),
               Line2D([0],[0],color="0.4",lw=1.6,ls="-",label="Tokyo transfer (no re-fit)"),
               Line2D([0],[0],color="0.4",lw=1.6,ls="--",label="Saitama re-estimation"),
               Line2D([0],[0],marker="o",color="none",markerfacecolor="0.5",markeredgecolor="none",markersize=5,label="Observed")]
    fig.legend(handles=handles, loc="lower center", ncol=5, fontsize=8.6, frameon=False, bbox_to_anchor=(0.5,-0.01))
    fig.suptitle(f"{POINT.capitalize()}: transfer of the legacy and revised models, and re-estimation of the revised model", fontsize=12.5)
    fig.tight_layout(rect=[0, 0.04, 1, 0.96])
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_path}")

def plot_seasonal(results, out_path):
    species = ["Ae", "Cx"]
    panels = [("B_PCMP", "#1f4e79", "-",  "B_PCMP\n(Tokyo transfer)"),
              ("MM",     "#c0392b", "-",  "MM\n(Tokyo transfer)"),
              ("MM_S",   "#c0392b", "--", "MM\n(Saitama re-estimation)")]
    fig, axes = plt.subplots(2, 3, figsize=(13.5, 7.2), sharex=True, sharey="row")
    wk_axis = np.arange(1, 53)
    for ri, sk in enumerate(species):
        info = SPECIES_INFO[sk]
        R = results[sk]; obs = R["obs"]; wl = len(obs)
        woy = np.array([(i % 52) + 1 for i in range(wl)])
        ob_mean = np.full(52, np.nan); ob_sd = np.full(52, np.nan)
        for w in range(52):
            sel = (woy == (w + 1)) & (obs >= 0.0)
            if np.any(sel):
                vals = np.log10(obs[sel] + 1.0); ob_mean[w] = np.mean(vals); ob_sd[w] = np.std(vals)
        def wkmean(mk):
            pm = np.full(52, np.nan)
            est = R["est"][mk]
            for w in range(52):
                sel = (woy == (w + 1))
                if np.any(sel): pm[w] = np.mean(np.log10(est[sel] + 1.0))
            return pm
        for ci, (mk, color, ls, title) in enumerate(panels):
            ax = axes[ri, ci]
            ax.fill_between(wk_axis, ob_mean - ob_sd, ob_mean + ob_sd, color="0.8", alpha=0.55, lw=0, zorder=1)
            ax.scatter(wk_axis, ob_mean, s=12, color="0.4", zorder=2, linewidths=0)
            if mk in R["est"]:
                ax.plot(wk_axis, wkmean(mk), color=color, ls=ls, lw=1.8, zorder=3)
            ax.grid(True, ls=":", alpha=0.4); ax.set_axisbelow(True); ax.set_xlim(1,52)
            if ri == 0: ax.set_title(title, fontsize=10.5)
            if ci == 0:
                g, s = info['pretty'].split()[0], info['pretty'].split()[1]
                ax.set_ylabel(f"$\\it{{{g}}}$ $\\it{{{s}}}$\nmean log10 (N + 1)", fontsize=9.5)
            if ri == 1: ax.set_xlabel("Week of year", fontsize=9.5)
    from matplotlib.lines import Line2D
    handles = [Line2D([0],[0],color="#1f4e79",lw=2,label="B_PCMP (legacy)"),
               Line2D([0],[0],color="#c0392b",lw=2,label="MM (revised, best model)"),
               Line2D([0],[0],color="0.4",lw=1.8,ls="-",label="Tokyo transfer (no re-fit)"),
               Line2D([0],[0],color="0.4",lw=1.8,ls="--",label="Saitama re-estimation"),
               Line2D([0],[0],marker="o",color="none",markerfacecolor="0.4",markeredgecolor="none",markersize=5,label="Observed mean")]
    fig.legend(handles=handles, loc="lower center", ncol=5, fontsize=8.6, frameon=False, bbox_to_anchor=(0.5,-0.01))
    fig.suptitle(f"{POINT.capitalize()} seasonal mean: transfer of the legacy and revised models, and re-estimation of the revised model", fontsize=12.5)
    fig.tight_layout(rect=[0, 0.04, 1, 0.96])
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_path}")

def main():
    global POINT, OBSPOINT
    ap = argparse.ArgumentParser()
    ap.add_argument("--param-dir", default=PARAM_DIR)
    ap.add_argument("--data-base", default=DATA_BASE)
    ap.add_argument("--point", default=POINT)
    ap.add_argument("--obs-point", default=OBSPOINT)
    ap.add_argument("--start-year", type=int, default=START_YEAR)
    ap.add_argument("--year-length", type=int, default=YEAR_LENGTH)
    ap.add_argument("--warmup-start", type=int, default=WARMUP_START)
    ap.add_argument("--warmup-end", type=int, default=WARMUP_END)
    ap.add_argument("--integrator", default=INTEGRATOR, choices=["E", "R"])
    
    # CLI overrides for parameter JSONs (B_PCMP and MM; MM_S = Saitama re-estimation)
    ap.add_argument("--ae-bpcmp", default=None)
    ap.add_argument("--ae-mm", default=None)
    ap.add_argument("--ae-mm-s", default=None)
    ap.add_argument("--cx-bpcmp", default=None)
    ap.add_argument("--cx-mm", default=None)
    ap.add_argument("--cx-mm-s", default=None)
    args = ap.parse_args()

    POINT = args.point
    OBSPOINT = args.obs_point

    explicit = {
        ("Ae", "B_PCMP"): args.ae_bpcmp, ("Ae", "MM"): args.ae_mm,
        ("Ae", "MM_S"): args.ae_mm_s,
        ("Cx", "B_PCMP"): args.cx_bpcmp, ("Cx", "MM"): args.cx_mm,
        ("Cx", "MM_S"): args.cx_mm_s,
    }

    results = {}
    metric_rows = []
    for sk in ("Ae", "Cx"):
        info = SPECIES_INFO[sk]
        print(f"\n=== {info['pretty']} : loading Saitama data ===")
        arrays = load_input_arrays(
            info["code"], info["full"], info["obs_candidates"], args.data_base,
            args.point, args.obs_point, args.start_year, args.year_length,
            args.warmup_start, args.warmup_end)
        xaxis = build_week_year_axis(arrays["data_year"], args.start_year, arrays["week_length"])
        
        R = {"obs": arrays["obs"], "xaxis": xaxis, "est": {}, "metrics": {}}
        
        for mk_full in MODEL_ORDER:
            is_saitama = mk_full.endswith("_S")
            base_mk = mk_full.replace("_S", "")
            target_region = "saitama" if is_saitama else "tokyo"
            
            path = explicit.get((sk, mk_full)) or PARAM_FILES.get((sk, mk_full))
            if path is None:
                path = autodiscover(args.param_dir, info["code"], base_mk, args.integrator, target_region)
                
            if path is None or not os.path.exists(path):
                print(f"  [{sk} {mk_full}] parameter JSON not found — skipping.")
                continue
                
            x, wash_code, use_rk4, label = x_from_json(path)
            print(f"  [{sk} {mk_full}] using {os.path.basename(path)} "
                  f"(washout={[k for k,v in WASH_CODE.items() if v==wash_code][0]}, "
                  f"{'RK4' if use_rk4 else 'Euler'})")
                  
            est = simulate_est_week(arrays, x, info["is_ae"], use_rk4,
                                    enable_diapause=1, enable_ccl=1, wash_code=wash_code)
            R["est"][mk_full] = est
            mt = metrics(arrays["obs"], est)
            R["metrics"][mk_full] = mt
            metric_rows.append({"species": info["pretty"], "model": mk_full,
                                "label": label, **mt})
            print(f"      -> R2_log={mt['R2_Log10']:.3f}  RMSE_log={mt['RMSE_Log10']:.3f}")
            
        results[sk] = R

    # 図
    plot_timeseries(results, args.start_year, f"fig5_transfer_reestimation_timeseries_{POINT}.pdf")
    plot_seasonal(results, f"fig5_transfer_reestimation_seasonal_{POINT}.pdf")

    # メトリクス CSV
    if metric_rows:
        dfm = pd.DataFrame(metric_rows)
        csv_path = f"transfer_vs_recalc_metrics_{POINT}.csv"
        dfm.to_csv(csv_path, index=False)
        print(f"\nSaved: {csv_path}")
        print("\n=== Metrics Summary ===")
        print(dfm.to_string(index=False))
    else:
        print("\nNo metrics computed. Check JSON files in directory.")

if __name__ == "__main__":
    main()