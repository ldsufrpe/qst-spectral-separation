"""Figure sweep (WORKFLOW Tarefa 4 `itf_vs_bounds`, finding R-B5), K4 family.

Panel a: Delta_mod = 24, k = 3..12 (D = 2k-1): exact ITF S^2 vs Sigma1^2, b_eta^2, B_E^2.
Panel b: k = 5, Delta_mod chosen so that g/(2 delta) takes the values in RATIOS,
         down to the edge g -> 2 delta^+ where the method stops informing.
kappa is optimized on the grid 1/1000 for Sigma1 (illustrative, not a table value).
"""
import csv
import sys
import time
from fractions import Fraction as Q

import mpmath as mp

import bounds as BD
import logutil as LU
from families import build
from run_benchmark import adaptive_instance, fmt, PROV

RES = LU.EXP / "results"
RATIOS = [Q(100), Q(30), Q(10), Q(5), Q(3), Q(2), Q(3, 2), Q(6, 5), Q(21, 20)]


def cliques(k, Delta):
    G = build("cliques", k)
    G["h"] = [Delta * G["comm"][v] for v in range(G["n"])]
    return G


def point(panel, G, Delta, ratio=None):
    t0, ts = time.time(), LU.now()
    mp.mp.dps = 50
    par = BD.derived(G, BD.block_data(G))
    kap = BD.best_kappa(G, par, den=1000)
    kq = Q(int(mp.nint(kap * 1000)), 1000)
    out, checks, _ = adaptive_instance(G, kq, 60)
    row = dict(panel=panel, k=G["k"], n=G["n"], D=out["D"], Delta=float(Delta),
               g_over_2delta=float(ratio) if ratio is not None else float((Delta - 4) / Q(1, 5)),
               eta=fmt(par["eta"]), kappa=str(kq), gamma_kappa=fmt(out["gamma_kappa"]),
               ITF=fmt(out["ITF"]), Fbar=fmt(out["Fbar"]), F_samp=fmt(out["F_samp"]),
               Sigma1_sq=fmt(min(1, out["Sigma1"] ** 2)), b_eta2=fmt(out["b_eta2"]),
               B_E2=fmt(out["B_E2"]), B2=fmt(out["B2"]),
               ratio_B2_over_ITF=fmt(out["ratio_B2_over_ITF"]), dps=out["dps"],
               checks_all_pass=all(checks.values()))
    print(" ".join(f"{k}={v}" for k, v in row.items()), flush=True)
    return row, checks, ts, t0


def main():
    LU.capture_env()
    hdr = LU.header("2026-09-18_exact_benchmarks", "numerical", ["numerical"], PROV,
                    "2026-09-18T00:00:00-03:00")
    jobs = [("a", cliques(k, Q(24)), Q(24), None) for k in range(3, 13)]
    jobs += [("b", cliques(5, 4 + r / 5), 4 + r / 5, r) for r in RATIOS]
    rows, bad = [], 0
    for panel, G, Delta, r in jobs:
        row, checks, ts, t0 = point(panel, G, Delta, r)
        rows.append(row)
        ok = all(checks.values())
        bad += not ok
        LU.append_run(hdr, dict(
            timestamp_start=ts, timestamp_end=LU.now(), wallclock_seconds=round(time.time() - t0, 1),
            peak_memory_mb=round(LU.peak_mb(), 1), status="success" if ok else "failure", seed=None,
            parameters=dict(sweep_panel=panel, family="cliques", k=G["k"], Delta_mod=str(Delta),
                            g_over_2delta=str(r) if r is not None else "", kappa=row["kappa"]),
            outputs=dict(scalar_results=row | {"checks": checks}),
            exit_code=0 if ok else 1, provenance_note=PROV,
            notes="seed=null: deterministic. Sweep point for figures/itf_vs_bounds.pdf."))
    with open(RES / "sweep.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
