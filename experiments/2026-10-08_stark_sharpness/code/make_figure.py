"""Figure: sharpness of the chain bound on the uniform Wannier-Stark ladder.

log10(Sigma/S) against N at F/J = 10: lines on Z for Sigma1, N1 and the three N2' forms,
with the asymptotic law log10(pi N/2) + N log10 B(q'^2) - log10((1-4q'^2)(1-eta^2)) of the
per-bridge form; finite chains (N2 and the rigorous N3x) as separate markers.
Single panel since Phase 3 S9 (2026-10-13, the author's choice): the former panel (b), fitted
against predicted c over F/J, was dropped; the law at every F/J is stated in the text with
its spread (results/fits.csv, N2p_bridge_S_level).
Writes figures/2026-10-08_stark_sharpness_rate.{pdf,png}, figures/figure_manifest.json,
results/viz_schema.json. IOP class not installed, so the J. Phys. A text width was not measured.
"""
import csv
import json
import math
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import scienceplots  # noqa: F401

EXP = Path(__file__).resolve().parent.parent
EID = "2026-10-08_stark_sharpness"
DESC = "rate"
FJ = 10                                              # F/J of the plotted ladder (J = 1)
# LEGEND LABELS -- the only place to edit. Set in Phase 3 S9 (2026-10-13) to the names the
# manuscript gives these results (Section 5 and Example ex:stark). figure-reviewer (Phase 5.6)
# checks legend against text.
LABELS = {
    "Sigma1": "Combes--Thomas bound",
    "N1": "bridge bound",
    "N2p_stated": r"chain bound (i), $\delta_{jl}=2\mu$",
    "N2p_ladder": r"chain bound (ii), $\delta_{jl}=2\mu$",
    "N2p_bridge": r"chain bound (ii), $\delta_{jl}=\mu$",
    "N2_chain": "walk-sum bound, finite chain",
    "N3x_chain": "contour certificate, finite chain",
    "pred": "asymptotic law",
}
ZSTYLE = [("Sigma1", "#E69F00", "s", "-"), ("N1", "#56B4E9", "^", "--"),
          ("N2p_stated", "#009E73", "D", "-."), ("N2p_ladder", "#0072B2", "P", ":"),
          ("N2p_bridge", "#D55E00", None, "-")]
ZSER = [(c, LABELS[c], col, mk, ls) for c, col, mk, ls in ZSTYLE]


def B(z):
    return (1 - math.sqrt(1 - 4 * z)) / (2 * z)


def main():
    z = [r for r in csv.DictReader(open(EXP / "results" / "z_series.csv")) if r["F_over_J"] == str(FJ)]
    ch = [r for r in csv.DictReader(open(EXP / "results" / "chains.csv")) if r["F_over_J"] == str(FJ)]
    plt.style.use(["science", "ieee"])
    matplotlib.rcParams.update({"figure.figsize": (4.22, 2.8), "savefig.dpi": 300,
                                "savefig.bbox": None,
                                "font.size": 8, "axes.labelsize": 8, "legend.fontsize": 7,
                                "xtick.labelsize": 7, "ytick.labelsize": 7,
                                "lines.linewidth": 1.0, "lines.markersize": 3.5})
    fig, ax = plt.subplots(constrained_layout=True)
    N = [int(r["N"]) for r in z]
    for c, lab, col, mk, ls in ZSER:
        ys = [float(r[f"{c}_log10_Sigma"]) - float(r["log10_S"]) for r in z]
        ax.plot(N, ys, color=col, linestyle=ls, marker=mk, markevery=(10, 10), label=lab)
    qp = 2 / FJ                                      # q' = 2J/F
    eta = 2 / (FJ - 2)                               # eta = delta/(g - delta), delta = 2J, g = F
    const = -math.log10((1 - 4 * qp ** 2) * (1 - eta ** 2))
    cp = math.log10(B(qp ** 2))
    ax.plot(N, [math.log10(math.pi * n / 2) + n * cp + const for n in N], color="#000000",
            linestyle=(0, (4, 2)), linewidth=0.8, label=LABELS["pred"])
    chs = sorted(ch, key=lambda r: int(r["N"]))
    ax.plot([int(r["N"]) for r in chs], [float(r["N2_log10_Sigma"]) - float(r["log10_S"]) for r in chs],
            linestyle="none", marker="o", color="#000000", markerfacecolor="none",
            label=LABELS["N2_chain"])
    c3 = [r for r in chs if r.get("N3x_log10_Sigma")]
    ax.plot([int(r["N"]) for r in c3], [float(r["N3x_log10_Sigma"]) - float(r["log10_S"]) for r in c3],
            linestyle="none", marker="X", color="#CC79A7", label=LABELS["N3x_chain"])
    ax.set_xlabel(r"number of bridges $N$ between $s$ and $t$")
    ax.set_ylabel(r"$\log_{10}(\Sigma_1/S)$")
    ax.set_ylim(-0.5, 32)
    ax.legend(loc="center right", frameon=False)
    out = EXP / "figures" / f"{EID}_{DESC}"
    fig.savefig(out.with_suffix(".pdf"))
    fig.savefig(out.with_suffix(".png"), dpi=200)
    manifest = {DESC: {
        "panels": [f"single panel: log10(Sigma/S) versus N at F/J = {FJ}"],
        "elements": [
            *[{"artist": "line" + ("+marker every 10 N from N = 11" if mk else ""),
               "encodes": f"{c} certificate on the infinite ladder (Z), log10(Sigma/S)",
               "visual": f"{col}, linestyle {ls!r}, marker {mk!r}", "panel": "all"}
              for c, _, col, mk, ls in ZSER],
            {"artist": "dashed line",
             "encodes": "asymptotic law log10(pi N/2) + N log10 B(q'^2) - log10((1-4q'^2)(1-eta^2)), "
                        "q' = 2J/F, eta = 2J/(F-2J), of the per-bridge form",
             "visual": "black, dashed (4,2), thin", "panel": "all"},
            {"artist": "hollow circles", "encodes": "N2 (Prop. N2, delta_jl = J) on finite chains L = N+1, against the chain's exact S",
             "visual": "black hollow circle, no line", "panel": "all"},
            {"artist": "X markers", "encodes": "rigorous N3x upper sum (M = 2000) on finite chains L <= 21",
             "visual": "#CC79A7 X, no line", "panel": "all"},
            {"artist": "legend", "encodes": "the eight series", "visual": "center right, no frame", "panel": "all"}],
        "axes": {"x": "N, the number of bridges between s and t", "y": "log10(Sigma/S)",
                 "xscale": "linear", "yscale": "linear (values are log10)"},
        "params": {"J": 1, "F_over_J": FJ, "N_range": "1..100 on Z; chains L = 6, 11, 21, 41, 61",
                   "figsize_in": [4.22, 2.8],
                   "legend_labels": "manuscript names, set in Phase 3 S9 (2026-10-13) via LABELS in code/make_figure.py",
                   "width_note": "iopjournal text width 153 mm = 6.02 in, measured in S11 (2026-10-14); printed at 0.7\\linewidth, no tight bbox"},
        "reference_lines": ["asymptotic law, dashed"],
        "normalization": "Sigma divided by the exact S of the same instance (Z: Bessel; chains: diagonalization)"}}
    (EXP / "figures" / "figure_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    viz = {"figures": [
        {"title": "Distance to exact S, log10(Sigma/S), vs N on Z (F/J = 10)", "plot_type": "line",
         "x_axis": "N", "y_axis": "dist_S", "z_axis": None, "series": "certificate",
         "filters": {"part": "Z", "F_over_J": "10"}, "log_scale_x": False, "log_scale_y": False},
        {"title": "Per-bridge N2' distance to S vs N on Z, by F/J", "plot_type": "line",
         "x_axis": "N", "y_axis": "dist_S", "z_axis": None, "series": "F_over_J",
         "filters": {"part": "Z", "certificate": "N2p_bridge"}, "log_scale_x": False, "log_scale_y": False},
        {"title": "Finite chains: distance to exact S vs N (F/J = 10)", "plot_type": "line",
         "x_axis": "N", "y_axis": "dist_S", "z_axis": None, "series": "certificate",
         "filters": {"part": "chain", "F_over_J": "10"}, "log_scale_x": False, "log_scale_y": False}]}
    (EXP / "results" / "viz_schema.json").write_text(json.dumps(viz, indent=2) + "\n")
    print("wrote", out.with_suffix(".pdf"))


if __name__ == "__main__":
    main()
