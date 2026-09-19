"""Variant for R-B4: fixed pair (s, t), growing tree.  This is the regime in which the
k-independence of the Theorem C prefactor matters: D stays fixed while k -> infinity.

s = vertex 3 of the root; t = vertex 3 of module 255 (leftmost node at depth 8), D = 17,
where the Theorem C bound is already informative.  Depth r = 8..16 (k up to 131071).
Certificates only: n >= 2044, beyond the multiprecision budget.
Uses the (kappa, beta, a, calC) chosen by run.py (results/theoremC_choice.csv).
"""
import csv
import math
import sys
import time
from fractions import Fraction as Q

import benchlib
import logutil as LU
from families import bfs
from run import consts, sigma1_float, RHO, ETA, M0, D0, PROV
from tree import build_tree

RES = LU.EXP / "results"


def main():
    ch = next(csv.DictReader(open(RES / "theoremC_choice.csv")))
    beta, calC = float(ch["beta"]), float(ch["calC"])
    LU.capture_env()
    hdr = LU.header("2026-09-18_binary_tree_family", "numerical", ["numerical", "networks"],
                    PROV, "2026-09-18T00:00:00-03:00")
    rows = []
    for r in range(8, 17):
        t0, ts = time.time(), LU.now()
        G = build_tree(r)
        G["t"] = 4 * 255 + 3
        ds_all, dt_all = bfs(G["adj"], G["s"]), bfs(G["adj"], G["t"])
        ds = [min(ds_all[v] for v in V) for V in G["parts"]]
        dt = [min(dt_all[v] for v in V) for V in G["parts"]]
        D = ds_all[G["t"]]
        assert D == 17
        kb = min((j / 1000 for j in range(1, 1000)), key=lambda x: sigma1_float(G, x, ds, dt))
        S1 = sigma1_float(G, kb, ds, dt)
        g2, c2, _ = consts(kb)
        C1u = ((2 * c2 * RHO * math.sqrt(M0) / math.sqrt(1 - ETA ** 2)
                + c2 ** 2 * RHO ** 2 * G["k"] * M0 / (1 - ETA ** 2)) * math.exp(g2 * D0))
        row = dict(r=r, k=G["k"], n=G["n"], D=D, kappa_sigma1=kb, Sigma1_sq=min(1.0, S1 ** 2),
                   C1unif=C1u, C1unif_bound_sq=min(1.0, (C1u * math.exp(-g2 * D)) ** 2),
                   C1unif_bound_sq_unclipped=(C1u * math.exp(-g2 * D)) ** 2,
                   thmC_bound=min(1.0, calC * math.exp(-2 * beta * D)), ITF="")
        ok = True
        rows.append(row | {"checks_all_pass": ok})
        print(row, flush=True)
        LU.append_run(hdr, dict(
            timestamp_start=ts, timestamp_end=LU.now(), wallclock_seconds=round(time.time() - t0, 1),
            peak_memory_mb=round(LU.peak_mb(), 1), status="success" if ok else "failure", seed=None,
            parameters=dict(phase="fixed_pair", r=r, t_module=255, beta=beta, kappa_sigma1=kb),
            outputs=dict(scalar_results={k_: (v if isinstance(v, (int, float, str)) else str(v))
                                         for k_, v in row.items()}),
            exit_code=0 if ok else 1, provenance_note=PROV,
            notes="seed=null: deterministic. Fixed (s,t), D=17, growing k."))
    with open(RES / "fixed_pair.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    sys.exit(main())
