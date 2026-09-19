"""Decay-rate analysis of the exact benchmark (R-B5): does the certified exponential
rate track the true decay of S with distance?

Two models for the exact S along each chain family (k modules, k-1 bridges):
  exponential   ln S = c0 - gamma_lin * D
  Stark ladder  ln S = c0 + (k-1) ln q - ln((k-1)!)
The second is the leading perturbative form when the site energies grow linearly
along the chain (energy denominators Delta, 2 Delta, 3 Delta, ...), which the
packing constraint (Prop. largura) makes unavoidable for (H2) with many modules.
Both are two-parameter least-squares fits; RMS residuals are compared.
Writes results/decay_fits.csv.
"""
import csv
import math
from pathlib import Path

import numpy as np

RES = Path(__file__).resolve().parent.parent / "results"


def fits(ks, Ds, lnS):
    ks, Ds, lnS = map(np.asarray, (ks, Ds, lnS))
    A = np.vstack([np.ones_like(Ds, float), -Ds]).T
    p_lin, *_ = np.linalg.lstsq(A, lnS, rcond=None)
    r_lin = lnS - A @ p_lin
    lf = np.array([math.lgamma(k) for k in ks])          # ln((k-1)!)
    B = np.vstack([np.ones_like(ks, float), ks - 1.0]).T
    p_st, *_ = np.linalg.lstsq(B, lnS + lf, rcond=None)
    r_st = lnS + lf - B @ p_st
    return dict(gamma_lin=p_lin[1], rms_lin=float(np.sqrt(np.mean(r_lin ** 2))),
                q_stark=float(np.exp(p_st[1])), rms_stark=float(np.sqrt(np.mean(r_st ** 2))),
                max_abs_resid_lin=float(np.max(np.abs(r_lin))),
                max_abs_resid_stark=float(np.max(np.abs(r_st))))


def main():
    rows = list(csv.DictReader(open(RES / "aggregate.csv")))
    series = {}
    for r in rows:
        series.setdefault(r["family"], []).append(
            (int(r["k"]), int(r["D"]), math.log(float(r["S"])) if float(r["S"]) > 0 else None,
             float(r["gamma_kappa"])))
    sweep = RES / "sweep.csv"
    if sweep.exists():
        pa = [r for r in csv.DictReader(open(sweep)) if r["panel"] == "a"]
        series["cliques_sweep_a"] = [(int(r["k"]), int(r["D"]), 0.5 * math.log(float(r["ITF"])),
                                      float(r["gamma_kappa"])) for r in pa]
    out = []
    for fam, pts in series.items():
        pts = [p for p in pts if p[2] is not None]
        if len(pts) < 3:
            continue
        f = fits([p[0] for p in pts], [p[1] for p in pts], [p[2] for p in pts])
        out.append(dict(series=fam, n_points=len(pts), k_range=f"{pts[0][0]}-{pts[-1][0]}",
                        gamma_true_linear_fit=round(f["gamma_lin"], 4),
                        gamma_kappa_certified_range=f"{min(p[3] for p in pts):.4f}-{max(p[3] for p in pts):.4f}",
                        rms_resid_exponential=round(f["rms_lin"], 4),
                        rms_resid_stark=round(f["rms_stark"], 4),
                        max_resid_exponential=round(f["max_abs_resid_lin"], 4),
                        max_resid_stark=round(f["max_abs_resid_stark"], 4),
                        q_stark=f"{f['q_stark']:.4e}"))
        print(out[-1])
    with open(RES / "decay_fits.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0].keys()))
        w.writeheader()
        w.writerows(out)


if __name__ == "__main__":
    main()
