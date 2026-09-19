"""figures/itf_vs_bounds.pdf  (WORKFLOW Tarefa 4; finding R-B5).

Reads results/sweep.csv.  Panel (a): K4 family, Delta_mod = 24, versus D.
Panel (b): k = 5, versus g/(2 delta) down to 1^+.  Axis label says ITF, never F_max.
Also writes figures/figure_manifest.json (every artist drawn below).
"""
import csv
import json
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import scienceplots  # noqa: F401
from cycler import cycler

EXP = Path(__file__).resolve().parent.parent
plt.style.use(["science", "ieee"])
matplotlib.rcParams.update({
    "figure.figsize": (6.48, 2.6),          # elsevier-1col full width (LAA), two panels
    "savefig.dpi": 300, "savefig.bbox": "tight", "savefig.pad_inches": 0.05,
    "savefig.format": "pdf",
    "axes.prop_cycle": cycler(color=["#000000", "#E69F00", "#56B4E9", "#009E73",
                                     "#F0E442", "#0072B2", "#D55E00", "#CC79A7"]),
})
OI = ["#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"]

SERIES = [  # (column, label, colour, marker, linestyle)
    ("ITF", r"ITF $S^2$ (exact)", OI[0], "o", "-"),
    ("Sigma1_sq", r"$\min\{1,\Sigma_1^2\}$", OI[1], "s", "-"),
    ("b_eta2", r"$b_\eta^2$", OI[2], "^", "--"),
    ("B_E2", r"$B_E^2$", OI[3], "D", "-."),
]


def main():
    rows = list(csv.DictReader(open(EXP / "results" / "sweep.csv")))
    pa = sorted((r for r in rows if r["panel"] == "a"), key=lambda r: int(r["D"]))
    pb = sorted((r for r in rows if r["panel"] == "b"), key=lambda r: float(r["g_over_2delta"]))
    fig, (axa, axb) = plt.subplots(1, 2)
    for col, lab, colr, mk, ls in SERIES:
        axa.plot([int(r["D"]) for r in pa], [float(r[col]) for r in pa], color=colr, marker=mk,
                 ls=ls, ms=3, label=lab)
        axb.plot([float(r["g_over_2delta"]) for r in pb], [float(r[col]) for r in pb], color=colr,
                 marker=mk, ls=ls, ms=3, label=lab)
    axa.set_yscale("log")
    axa.set_xlabel(r"distance $D=\mathrm{dist}(s,t)$")
    axa.set_ylabel(r"ITF and certificates")
    axa.set_title(r"(a) $\Delta_{\rm mod}=24$, $k=3,\dots,12$", fontsize=8)
    axb.set_xscale("log")
    axb.set_yscale("log")
    axb.axvline(1.0, color="0.5", ls=":", lw=0.8)
    axb.set_xlabel(r"$g/(2\delta)$")
    axb.set_title(r"(b) $k=5$, $D=9$", fontsize=8)
    axa.legend(loc="lower left", fontsize=6, frameon=False)
    fig.tight_layout()
    out = EXP / "figures" / "itf_vs_bounds.pdf"
    fig.savefig(out)
    fig.savefig(EXP / "figures" / "itf_vs_bounds.png", format="png", dpi=200)

    manifest = {"itf_vs_bounds": {
        "panels": ["(a) K4-module chain, Delta_mod = 24, k = 3..12, versus D",
                   "(b) K4-module chain, k = 5 (D = 9), versus g/(2 delta)"],
        "elements": [
            {"artist": "line+circle markers", "encodes": "exact ITF S(s,t)^2 (multiprecision diagonalization)",
             "visual": "black, solid, circle", "panel": "all"},
            {"artist": "line+square markers", "encodes": "spatial certificate min{1, Sigma1^2} at grid-optimal kappa",
             "visual": "orange #E69F00, solid, square", "panel": "all"},
            {"artist": "line+triangle markers", "encodes": "subspace bound b_eta^2 (Cor. comparacao)",
             "visual": "sky blue #56B4E9, dashed, triangle", "panel": "all"},
            {"artist": "line+diamond markers", "encodes": "energy bound B_E^2 (Prop. energia)",
             "visual": "bluish green #009E73, dash-dot, diamond", "panel": "all"},
            {"artist": "legend", "encodes": "the four series", "visual": "lower left of (a), no frame", "panel": "(a)"},
            {"artist": "panel titles", "encodes": "parameters of each panel", "visual": "8 pt text above axes", "panel": "all"},
        ],
        "axes": {"x": "(a) graph distance D; (b) g/(2 delta), dimensionless",
                 "y": "ITF and certificate values (dimensionless, shared meaning, not shared axis)",
                 "xscale": "(a) linear; (b) log", "yscale": "log"},
        "params": {"family": "K4 modules (Example 7.1)", "w": 0.1, "m": 4,
                   "panel_a": "Delta_mod = 24, k = 3..12",
                   "panel_b": "k = 5, Delta_mod = 4 + 0.2 * (g/2delta), g/(2delta) in {1.05,1.2,1.5,2,3,5,10,30,100}",
                   "kappa": "grid-optimal for Sigma1 on step 1e-3 (1/2 on plateaus)",
                   "s,t": "first vertex of V_0, last vertex of V_{k-1}"},
        "reference_lines": ["(b) vertical dotted grey line at g/(2 delta) = 1, the edge of hypothesis (H2)"],
        "normalization": "none; all quantities are probabilities or bounds on them in [0,1]",
    }}
    (EXP / "figures" / "figure_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print("wrote", out)


if __name__ == "__main__":
    main()
