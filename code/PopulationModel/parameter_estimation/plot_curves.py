#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_curves.py
==============
Stage-mortality curve plotter (Aedes original vs Culex re-estimated), generated
directly from one or more *_params.json files written by main_cx.py /
main_cx_saitama.py.

This is a STANDALONE re-implementation of main_cx.py's plot_ae_vs_cx_curves():
it does NOT re-run estimation, it just reads the vertex thermal knobs stored in
the JSON and draws the three stage curves (larva, pupa, adult).

Two modes
---------
1) ONE json  -> a single 3-panel figure (Aedes vs that Culex fit), exactly like
   the figure main_cx.py emits, but reproducible from the saved JSON.

2) TWO jsons (e.g. Tokyo and Saitama) -> the same 3 panels, but with BOTH Culex
   fits overlaid on the Aedes baseline, so the site-to-site difference in the
   estimated thermal curves is visible at a glance.

Usage
-----
# single fit
python3 plot_curves.py --json Cx_MEM+Ta_saitama_params.json

# compare two fits (any number of --json up to a few; first is drawn solid,
# the rest dashed/dotted). Give each a label with --label in the same order.
python3 plot_curves.py \
    --json Cx_MEM+Ta_tokyo_params.json   --label "Tokyo (MM+Ta)" \
    --json Cx_MEM+Ta_saitama_params.json --label "Saitama (MM+Ta)" \
    --out  Cx_curves_tokyo_vs_saitama

# Saitama-only convenience (just a normal single-json call):
python3 plot_curves.py --json Cx_MEM+Ta_saitama_params.json --out Cx_curves_saitama
"""

import os
import json
import argparse
import numpy as np

# ---------------------------------------------------------------------------
# Constants copied verbatim from main_cx.py (Jia 2016 Aedes vertex parameters
# and the literature anchors). Keep these in sync with main_cx.py if it changes.
# ---------------------------------------------------------------------------
M1_TW_MAX = 35.04
M3_TA_MAX = 39.2
M1_TW_MIN = 5.16
M2_TW_MIN = 5.16
M3_TA_MIN = 3.18

AE_LARVA_TOPT = 14.8238
AE_LARVA_MMIN = 0.016805
AE_LARVA_K    = 0.1305
AE_PUPA_TOPT  = 16.8322
AE_PUPA_MMIN  = 0.021740
AE_PUPA_K     = 0.1502
AE_ADULT_TOPT = 21.2103
AE_ADULT_MMIN = 0.015793
AE_ADULT_K    = 0.1921

CX_MORT_CAP    = 0.8
CX_TOPT_LIT_LO = 20.0
CX_TOPT_LIT_HI = 28.0
CX_LETHAL_LIT  = 35.0

# Thermal-knob parameter names in the JSON "parameters" block.
CX_THERMAL_NAMES = [
    "p_larva_dTopt", "p_larva_fwarm", "p_larva_fcool",
    "p_pupa_dTopt",  "p_pupa_fwarm",  "p_pupa_fcool",
    "p_adult_dTopt", "p_adult_fwarm", "p_adult_fcool",
]
CX_THERMAL_NEUTRAL = [0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0]

# Per-stage definition: (display name, Ae T_opt, m_min, Ae k, lo, hi, knob base, temp type)
# knob base indexes into the 9-vector [larva(0..2), pupa(3..5), adult(6..8)].
STAGE_DEFS = (
    ("Larva (m_L)", AE_LARVA_TOPT, AE_LARVA_MMIN, AE_LARVA_K, M1_TW_MIN, M1_TW_MAX, 0, "water"),
    ("Pupa (m_P)",  AE_PUPA_TOPT,  AE_PUPA_MMIN,  AE_PUPA_K,  M2_TW_MIN, M1_TW_MAX, 3, "water"),
    ("Adult (m_A)", AE_ADULT_TOPT, AE_ADULT_MMIN, AE_ADULT_K, M3_TA_MIN, M3_TA_MAX, 6, "air"),
)

# Colour cycle for multiple Culex fits.
CX_COLORS = ["#2471A3", "#6c3483", "#117A65", "#B9770E", "#922B21"]
CX_LSTYLE = ["--", "-.", ":", (0, (3, 1, 1, 1)), (0, (5, 1))]


def vertex_curve(T, topt, mmin, k_cool, k_warm, lo, hi, cap=CX_MORT_CAP):
    """Vectorised vertex-form stage mortality, identical to main_cx.py's curve()."""
    T = np.asarray(T, float)
    k = np.where(T < topt, k_cool, k_warm)
    denom = (1.0 / mmin) - k * (T - topt) ** 2
    m = np.where(denom > 1e-12, 1.0 / denom, cap)
    m = np.clip(m, 0.0, cap)
    m = np.where((T > lo) & (T < hi), m, cap)
    return m


def load_thermal_knobs(path):
    """Return (knobs9, label, scope) from a *_params.json file. Missing thermal
    params fall back to neutral (Aedes) so legacy JSON still plots as Aedes."""
    with open(path, "r", encoding="utf-8") as f:
        j = json.load(f)
    params = j.get("parameters", {})
    if all(nm in params for nm in CX_THERMAL_NAMES):
        knobs = [float(params[nm]) for nm in CX_THERMAL_NAMES]
    else:
        knobs = list(CX_THERMAL_NEUTRAL)
    label = j.get("label", os.path.splitext(os.path.basename(path))[0])
    scope = _scope_of(knobs)
    return knobs, label, scope


def _scope_of(knobs):
    def moved(base):
        return (abs(knobs[base] - 0.0) > 1e-6 or
                abs(knobs[base + 1] - 1.0) > 1e-6 or
                abs(knobs[base + 2] - 1.0) > 1e-6)
    larva, pupa, adult = moved(0), moved(3), moved(6)
    if not (larva or pupa or adult):
        return "none"
    if adult and not (larva or pupa):
        return "adult"
    return "all"


def plot_curves(fits, out_base, title_extra=""):
    """fits: list of (knobs9, label). Draw the 3-panel Ae-vs-Cx figure."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        print(f"matplotlib unavailable: {e}")
        return None

    Tg = np.linspace(0, 40, 500)
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
    for ax, (name, ae_topt, mmin, ae_k, lo, hi, base, ttype) in zip(axes, STAGE_DEFS):
        # Aedes baseline (neutral knobs).
        m_ae = vertex_curve(Tg, ae_topt, mmin, ae_k, ae_k, lo, hi)
        ax.plot(Tg, m_ae, color="#C0392B", lw=2.2, label="Ae (Jia 2016)", zorder=3)

        # Each Culex fit.
        topt_strs = []
        for fi, (knobs, lab) in enumerate(fits):
            dTopt = knobs[base]; fwarm = knobs[base + 1]; fcool = knobs[base + 2]
            cx_topt = ae_topt + dTopt
            m_cx = vertex_curve(Tg, cx_topt, mmin, ae_k * fcool, ae_k * fwarm, lo, hi)
            color = CX_COLORS[fi % len(CX_COLORS)]
            ls = CX_LSTYLE[fi % len(CX_LSTYLE)] if len(fits) > 1 else "--"
            cx_label = (lab if len(fits) > 1 else "Cx (re-estimated)")
            ax.plot(Tg, m_cx, color=color, lw=2.2, ls=ls, label=cx_label, zorder=4)
            topt_strs.append(f"{cx_topt:.1f}")

        # Literature shading / lethal line.
        ax.axvspan(CX_TOPT_LIT_LO, CX_TOPT_LIT_HI, color="#2ecc71", alpha=0.10,
                   label="Cx optimum (lit. 20-28C)")
        ax.axvline(CX_LETHAL_LIT, color="gray", ls=":", lw=1.2,
                   label="Cx ~lethal (lit. ~35C)")

        topt_join = " / ".join(topt_strs)
        ax.set_title(f"{name} [{ttype} temp]\nAe Topt={ae_topt:.1f} -> Cx Topt={topt_join} C")
        ax.set_xlabel(f"{ttype} temperature (C)")
        ax.set_ylabel("daily mortality /day")
        ax.set_ylim(0, 0.5)
        ax.grid(alpha=0.25)
        if base == 0:
            ax.legend(fontsize=7, loc="upper left")

    if len(fits) > 1:
        sup = "Stage mortality: Aedes (Jia) vs Culex (re-estimated) -- " + \
              " vs ".join(lab for _, lab in fits)
    else:
        sup = f"Stage mortality: Aedes (Jia) vs Culex (re-estimated) -- {fits[0][1]}"
    if title_extra:
        sup += f"  {title_extra}"
    fig.suptitle(sup, fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.94])

    png = out_base + ".png"
    pdf = out_base + ".pdf"
    fig.savefig(png, dpi=130)
    fig.savefig(pdf)
    plt.close(fig)
    print(f"Saved: {png}")
    print(f"Saved: {pdf}")

    # Print a small text summary (T_opt, warm-lethal) for each fit.
    for knobs, lab in fits:
        print(f"\n[{lab}]  thermal scope = {_scope_of(knobs)}")
        for name, ae_topt, mmin, ae_k, lo, hi, base, ttype in STAGE_DEFS:
            dTopt = knobs[base]; fwarm = knobs[base + 1]; fcool = knobs[base + 2]
            cx_topt = ae_topt + dTopt
            # warm-lethal: where the warm-side curve hits the cap.
            Tw = np.linspace(cx_topt, hi, 400)
            mw = vertex_curve(Tw, cx_topt, mmin, ae_k * fcool, ae_k * fwarm, lo, hi)
            hit = np.where(mw >= CX_MORT_CAP - 1e-9)[0]
            leth = Tw[hit[0]] if hit.size else hi
            print(f"    {name:11s}: T_opt {ae_topt:.1f}->{cx_topt:.1f} C  "
                  f"k_warm x{fwarm:.2f}  k_cool x{fcool:.2f}  warm-lethal ~{leth:.1f} C")
    return png


def main():
    ap = argparse.ArgumentParser(description="Plot Ae vs Cx stage-mortality curves from params JSON(s).")
    ap.add_argument("--json", action="append", required=True,
                    help="Path to a *_params.json (repeatable to overlay several fits).")
    ap.add_argument("--label", action="append", default=None,
                    help="Label for each --json, in the same order. Optional.")
    ap.add_argument("--out", default=None,
                    help="Output file base name (without extension). "
                         "Default derives from the first JSON's label.")
    args = ap.parse_args()

    fits = []
    for i, path in enumerate(args.json):
        if not os.path.exists(path):
            raise FileNotFoundError(path)
        knobs, label, scope = load_thermal_knobs(path)
        if args.label and i < len(args.label):
            label = args.label[i]
        fits.append((knobs, label))

    if args.out:
        out_base = args.out
    else:
        safe = fits[0][1].replace("[", "_").replace("]", "").replace("/", "_").replace(" ", "_")
        out_base = f"{safe}_AE_vs_CX_curves"

    plot_curves(fits, out_base)


if __name__ == "__main__":
    main()
