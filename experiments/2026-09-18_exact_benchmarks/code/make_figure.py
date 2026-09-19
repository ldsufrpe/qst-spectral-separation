"""figures/itf_vs_bounds.pdf  (WORKFLOW Tarefa 4; finding R-B5).

Reads results/sweep.csv.  Panel (a): K4 family, Delta_mod = 24, versus D.
Panel (b): k = 5, versus g/(2 delta) down to 1^+.  Axis label says ITF, never F_max.
Panels (c), (d) (added 2026-09-19, replacing Table 5 of the manuscript; no new
computation): read ../2026-09-18_binary_tree_family/results/aggregate.csv (root to a
leaf at depth r, versus D = 2r+1) and fixed_pair.csv (fixed pair, D = 17, versus the
number of modules k).
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
    "figure.figsize": (6.48, 5.0),          # elsevier-1col full width (LAA), 2x2 panels
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
TREE = EXP.parent / "2026-09-18_binary_tree_family" / "results"
TREE_SERIES = [  # (column, label, colour, marker, linestyle)
    ("ITF", r"ITF $S^2$ (exact)", OI[0], "o", "-"),
    ("Sigma1_sq", r"$\min\{1,\Sigma_1^2\}$", OI[1], "s", "-"),
    ("C1unif_bound_sq", r"diameter form", OI[6], "v", ":"),
    ("thmC_bound", r"Theorem C", OI[5], "P", "-"),
]


def main():
    rows = list(csv.DictReader(open(EXP / "results" / "sweep.csv")))
    pa = sorted((r for r in rows if r["panel"] == "a"), key=lambda r: int(r["D"]))
    pb = sorted((r for r in rows if r["panel"] == "b"), key=lambda r: float(r["g_over_2delta"]))
    fig, ((axa, axb), (axc, axd)) = plt.subplots(2, 2)
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
    axb.set_yticks([10.0**e for e in range(0, -35, -5)])  # decades anchored at 10^0
    axb.set_ylim(top=30.0)  # data lie in [0,1]; headroom so markers at 1 are not clipped
    axb.axvline(1.0, color="0.5", ls=":", lw=0.8)
    axb.set_xlabel(r"$g/(2\delta)$")
    axb.set_title(r"(b) $k=5$, $D=9$", fontsize=8)
    axa.legend(loc="lower left", fontsize=8, frameon=False)  # prints at 7.5 pt

    tc = list(csv.DictReader(open(TREE / "aggregate.csv")))
    td = list(csv.DictReader(open(TREE / "fixed_pair.csv")))
    for col, lab, colr, mk, ls in TREE_SERIES:
        pts = [(int(r["D"]), min(1.0, float(r[col]))) for r in tc if r[col]]
        axc.plot([x for x, _ in pts], [y for _, y in pts], color=colr, marker=mk, ls=ls, ms=3, label=lab)
        if col != "ITF":
            pts = [(int(r["k"]), min(1.0, float(r[col]))) for r in td]
            axd.plot([x for x, _ in pts], [y for _, y in pts], color=colr, marker=mk, ls=ls, ms=3, label=lab)
    axc.set_yscale("log")
    axc.set_xlabel(r"distance $D=2r+1$")
    axc.set_ylabel(r"ITF and certificates")
    axc.set_title(r"(c) tree, root to a leaf at depth $r=1,\dots,12$", fontsize=8)
    axc.set_yticks([10.0**e for e in range(0, -60, -10)])  # decades anchored at 10^0
    axc.set_ylim(top=30.0)
    axc.legend(loc="lower right", fontsize=8, frameon=False)  # empty region right of the ITF
    axd.set_xscale("log")
    axd.set_yscale("log")
    axd.set_xlabel(r"number of modules $k$")
    axd.set_title(r"(d) tree, fixed pair at $D=17$, $r=8,\dots,16$", fontsize=8)
    fig.tight_layout()
    out = EXP / "figures" / "itf_vs_bounds.pdf"
    fig.savefig(out)
    fig.savefig(EXP / "figures" / "itf_vs_bounds.png", format="png", dpi=200)

    manifest = {"itf_vs_bounds": {
        "panels": ["(a) K4-module chain, Delta_mod = 24, k = 3..12, versus D",
                   "(b) K4-module chain, k = 5 (D = 9), versus g/(2 delta)",
                   "(c) binary tree of K4 modules (Example 7.5), s in the root, t in a leaf at depth r = 1..12, versus D = 2r+1",
                   "(d) same tree, fixed pair with D = 17, depth r = 8..16, versus the number of modules k"],
        "elements": [
            {"artist": "line+circle markers", "encodes": "exact ITF S(s,t)^2 (multiprecision diagonalization)",
             "visual": "black, solid, circle", "panel": "all"},
            {"artist": "line+square markers", "encodes": "spatial certificate min{1, Sigma1^2} at grid-optimal kappa",
             "visual": "orange #E69F00, solid, square", "panel": "all"},
            {"artist": "line+triangle markers", "encodes": "subspace bound b_eta^2 (Cor. comparacao)",
             "visual": "sky blue #56B4E9, dashed, triangle", "panel": "all"},
            {"artist": "line+diamond markers", "encodes": "energy bound B_E^2 (Prop. energia)",
             "visual": "bluish green #009E73, dash-dot, diamond", "panel": "all"},
            {"artist": "line+down-triangle markers", "encodes": "diameter form min{1,(C1unif e^{-gamma D})^2}",
             "visual": "vermillion #D55E00, dotted, down-triangle", "panel": "(c),(d)"},
            {"artist": "line+plus markers", "encodes": "Theorem C bound min{1, calC e^{-2 beta D}}",
             "visual": "blue #0072B2, solid, filled plus", "panel": "(c),(d)"},
            {"artist": "legend", "encodes": "the four chain series", "visual": "lower left of (a), no frame", "panel": "(a)"},
            {"artist": "legend", "encodes": "the four tree series", "visual": "lower right of (c), no frame", "panel": "(c)"},
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
