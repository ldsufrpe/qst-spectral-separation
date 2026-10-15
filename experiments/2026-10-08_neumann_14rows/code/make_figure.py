"""Figure: certificates on S against k for the three Section 7 families (from aggregate.csv).

Writes figures/2026-10-08_neumann_14rows_bounds_vs_k.{pdf,png}, figures/figure_manifest.json
and results/viz_schema.json. Style: figures_style.md (SciencePlots, Okabe-Ito).
Width: 6.0 in full text width. The IOP class (iopart/iopjournal) is not installed here,
so the J. Phys. A text width was NOT measured; measure before submission.
"""
import csv
import json
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import scienceplots  # noqa: F401

EXP = Path(__file__).resolve().parent.parent
EID = "2026-10-08_neumann_14rows"
DESC = "bounds_vs_k"
FAMS = [("cliques", "(a) $K_4$ modules"), ("dimer", "(b) dimer modules"), ("p4", "(c) $P_4$ modules")]
# LEGEND LABELS -- the only place to edit. Set in Phase 3 S9 (2026-10-13) to the names the
# manuscript gives these results in Section 5; figure-reviewer (Phase 5.6) checks that legend
# and text agree. The hybrid bound N3 left the paper in the triage v3 (S10c, 2026-10-14);
# its rows stay in aggregate.csv.
LABELS = {
    "S": r"exact $S$",
    "floor": r"floor $\sum_i\|P_ie_s\|\|P_ie_t\|$",
    "Sigma1": "Combes--Thomas bound",
    "N1": "bridge bound",
    "N2": "walk-sum bound",
    "N2p_path": "chain bound",
    "N3x": "contour certificate",
}
# (certificate, colour, marker, linestyle) -- Okabe-Ito order, yellow skipped
STYLE = [
    ("S", "#000000", "o", "-"),
    ("floor", "#000000", None, ":"),
    ("Sigma1", "#E69F00", "s", "-"),
    ("N1", "#56B4E9", "^", "--"),
    ("N2", "#009E73", "D", "-."),
    ("N2p_path", "#0072B2", "P", ":"),
    ("N3x", "#CC79A7", "X", "-"),
]
SERIES = [(c, LABELS[c], col, mk, ls) for c, col, mk, ls in STYLE]


def main():
    rows = list(csv.DictReader(open(EXP / "results" / "aggregate.csv")))
    plt.style.use(["science", "ieee"])
    matplotlib.rcParams.update({"figure.figsize": (6.02, 2.35), "savefig.dpi": 300,
                                "savefig.bbox": None,
                                "font.size": 8, "axes.labelsize": 8, "legend.fontsize": 7,
                                "xtick.labelsize": 7, "ytick.labelsize": 7,
                                "lines.linewidth": 1.0, "lines.markersize": 3.5})
    fig, axes = plt.subplots(1, 3, constrained_layout=True)
    for ax, (fam, title) in zip(axes, FAMS):
        fr = [r for r in rows if r["family"] == fam]
        ks = sorted({int(r["k"]) for r in fr})
        for cert, label, col, mk, ls in SERIES:
            if cert in ("S", "floor"):
                key = "log10_S" if cert == "S" else "log10_floor"
                ys = [float(next(r[key] for r in fr if int(r["k"]) == k)) for k in ks]
            else:
                ys = [float(next(r["log10_sigma"] for r in fr
                                 if int(r["k"]) == k and r["certificate"] == cert)) for k in ks]
            ax.plot(ks, ys, color=col, marker=mk, linestyle=ls, label=label,
                    markerfacecolor="none" if cert == "S" else col, zorder=3 if cert == "S" else 2)
        ax.set_title(title, fontsize=8)
        ax.set_xlabel(r"number of modules $k$")
        ax.set_xticks(ks if len(ks) <= 6 else ks[::2])
    axes[0].set_ylabel(r"$\log_{10}$ of the bound on $S(s,t)$")
    h, lab = axes[0].get_legend_handles_labels()
    fig.legend(h, lab, loc="outside lower center", ncol=4, frameon=False)
    out = EXP / "figures" / f"{EID}_{DESC}"
    fig.savefig(out.with_suffix(".pdf"))
    fig.savefig(out.with_suffix(".png"), dpi=200)
    manifest = {DESC: {
        "panels": [f"{t}: family {f}, the Table 3 rows, versus k" for f, t in FAMS],
        "elements": [{"artist": "line+marker" if mk else "line",
                      "encodes": {"S": "exact S(s,t) (multiprecision diagonalization, 2026-09-18_exact_benchmarks)",
                                  "floor": "sum_i ||P_i e_s|| ||P_i e_t|| with exact block masses (best any block-mass certificate can give)",
                                  "Sigma1": "paper's Combes-Thomas certificate Sigma1 (Arb, 2026-09-18_interval_certification)",
                                  "N1": "Theorem N1 certificate (Arb)",
                                  "N2": "Proposition N2 community-graph walk sum (Arb)",
                                  "N2p_path": "Corollary N2', first-arrival path product (Arb)",
                                  "N3x": "rigorous N3x upper Riemann sum, min over M in {1000, 2000} (Arb)"}[c],
                      "visual": f"{col}, linestyle '{ls}', marker {mk!r}"
                                + (", hollow" if c == "S" else ""),
                      "panel": "all"} for c, _, col, mk, ls in SERIES]
        + [{"artist": "legend", "encodes": f"the {len(SERIES)} series", "visual": "below the panels, 4 columns, no frame", "panel": "all"},
           {"artist": "panel titles", "encodes": "family", "visual": "(a)-(c) above each panel", "panel": "all"}],
        "axes": {"x": "number of modules k", "y": "log10 of the value (each certificate bounds S from above)",
                 "xscale": "linear", "yscale": "linear (values are log10)"},
        "params": {"rows": "K4 k=3,4,5,6,8,10; dimer k=3,5,7,9,12; P4 k=10,20,40",
                   "certificate_precision": "Arb 128 bits, upper endpoints",
                   "figsize_in": [6.02, 2.35],
                   "legend_labels": "manuscript names, set in Phase 3 S9 (2026-10-13) via LABELS in code/make_figure.py",
                   "width_note": "iopjournal text width 153 mm = 6.02 in, measured in S11 (2026-10-14); printed at \\linewidth, no tight bbox"},
        "reference_lines": [],
        "normalization": "none; log10 of S and of each bound",
    }}
    (EXP / "figures" / "figure_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    viz = {"figures": [
        {"title": f"log10 of each certificate on S vs k ({f})", "plot_type": "line",
         "x_axis": "k", "y_axis": "log10_sigma", "z_axis": None, "series": "certificate",
         "filters": {"family": f}, "log_scale_x": False, "log_scale_y": False} for f, _ in FAMS
    ] + [
        {"title": f"Distance to the exact S, log10(Sigma/S), vs k ({f})", "plot_type": "line",
         "x_axis": "k", "y_axis": "dist_S", "z_axis": None, "series": "certificate",
         "filters": {"family": f, "certificate": ["Sigma1", "N1", "N2", "N2p_path", "N3", "N3x", "best"]},
         "log_scale_x": False, "log_scale_y": False} for f, _ in FAMS
    ]}
    (EXP / "results" / "viz_schema.json").write_text(json.dumps(viz, indent=2) + "\n")
    print("wrote", out.with_suffix(".pdf"))


if __name__ == "__main__":
    main()
