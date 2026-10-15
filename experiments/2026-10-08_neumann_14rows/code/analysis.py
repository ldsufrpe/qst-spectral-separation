"""MR2: share of the exact decay rate in k captured by each certificate, per family.

share = slope_k(log10 Sigma) / slope_k(log10 S), reported three ways:
  endpoint  - secant between the smallest and largest k (opportunities.md definition)
  lsq       - least-squares slope over all rows of the family
  lastpair  - secant between the two largest k (asymptotic proxy)
Reads results/aggregate.csv only; writes results/rate_shares.csv.
"""
import csv
from pathlib import Path

import numpy as np

EXP = Path(__file__).resolve().parent.parent
# opportunities.md round 2, "What the round-2 candidate list survived" (endpoint secant)
EXPECTED = {"Sigma1": (12, 17, 5), "N1": (41, 40, 6), "N2": (57, 59, 30), "N3vN2": (72, 83, 74)}
FAMS = ("cliques", "dimer", "p4")


def main():
    rows = list(csv.DictReader(open(EXP / "results" / "aggregate.csv")))
    certs = list(dict.fromkeys(r["certificate"] for r in rows))
    out = []
    for fam in FAMS:
        fr = [r for r in rows if r["family"] == fam]
        ks = sorted({int(r["k"]) for r in fr})
        S = np.array([float(next(r["log10_S"] for r in fr if int(r["k"]) == k)) for k in ks])
        kk = np.array(ks, float)
        sS = dict(endpoint=(S[-1] - S[0]) / (kk[-1] - kk[0]), lsq=np.polyfit(kk, S, 1)[0],
                  lastpair=(S[-1] - S[-2]) / (kk[-1] - kk[-2]))
        for c in certs:
            vals = {int(r["k"]): float(r["log10_sigma"]) for r in fr if r["certificate"] == c}
            if sorted(vals) != ks:
                continue
            y = np.array([vals[k] for k in ks])
            sb = dict(endpoint=(y[-1] - y[0]) / (kk[-1] - kk[0]), lsq=np.polyfit(kk, y, 1)[0],
                      lastpair=(y[-1] - y[-2]) / (kk[-1] - kk[-2]))
            rec = dict(family=fam, certificate=c, k_range=f"{ks[0]}-{ks[-1]}", n_points=len(ks))
            for d in ("endpoint", "lsq", "lastpair"):
                rec[f"slope_S_{d}"] = f"{sS[d]:.6f}"
                rec[f"slope_cert_{d}"] = f"{sb[d]:.6f}"
                rec[f"share_{d}"] = f"{sb[d] / sS[d]:.4f}"
            if c in EXPECTED:
                rec["expected_endpoint_pct"] = EXPECTED[c][FAMS.index(fam)]
            out.append(rec)
    with open(EXP / "results" / "rate_shares.csv", "w", newline="") as f:
        keys = []
        for r in out:
            keys += [x for x in r if x not in keys]
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(out)
    print(f"{'certificate':12s} " + "  ".join(f"{fam:>22s}" for fam in FAMS))
    print(f"{'':12s} " + "  ".join(f"{'end / lsq / last':>22s}" for _ in FAMS))
    for c in certs:
        cells = []
        for fam in FAMS:
            r = next((r for r in out if r["family"] == fam and r["certificate"] == c), None)
            cells.append("—".rjust(22) if r is None else
                         f"{100*float(r['share_endpoint']):5.1f} {100*float(r['share_lsq']):5.1f} "
                         f"{100*float(r['share_lastpair']):5.1f}".rjust(22))
        exp = EXPECTED.get(c)
        print(f"{c:12s} " + "  ".join(cells) + (f"   expected end {exp}" if exp else ""))


if __name__ == "__main__":
    main()
