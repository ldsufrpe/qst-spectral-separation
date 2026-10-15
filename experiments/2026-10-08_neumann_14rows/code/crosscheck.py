"""Verification step 3: the reduced N3x (bridge endpoints only) equals the dense N3x of
the toy, block by block, at the same M. Both evaluate the same exact quantity from the
same upper endpoints U, so they must agree to rounding (< 1e-25 relative).
Instances: three Table 3 rows (bridges form a matching) and a chain of singletons with
site energies 10n and couplings 1 (every interior vertex carries two bridges).

Usage: python crosscheck.py [--M 300]
"""
import argparse
import csv
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import flint                                         # noqa: E402
from flint import arb                                # noqa: E402
import neumann as nm                                 # noqa: E402
from flint import fmpq                               # noqa: E402


def singleton_chain(L, F=10, J=1):
    """Non-matching instance in the rig_setup format (as in the Stark experiment)."""
    F, J = fmpq(F), fmpq(J)
    absV = [[J if abs(u - v) == 1 else fmpq(0) for v in range(L)] for u in range(L)]
    return dict(name="singleton_chain", k=L, n=L, parts=[[v] for v in range(L)],
                comm=list(range(L)), lo=[v * F for v in range(L)], hi=[v * F for v in range(L)],
                W=[fmpq(0)] * L, gij=[[abs(i - j) * F if i != j else None for j in range(L)]
                                      for i in range(L)],
                g=F, delta=2 * J, rho=F / 2, eta=2 * J / (F - 2 * J), absV=absV,
                Vblk=[r[:] for r in absV], Ncomm=[[abs(i - j) for j in range(L)] for i in range(L)],
                Hblocks=[[[v * F]] for v in range(L)], B=list(range(L)), s=0, t=L - 1)

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--M", type=int, default=300)
    ap.add_argument("--prec", type=int, default=128)
    a = ap.parse_args()
    flint.ctx.prec = a.prec
    rows, worst = [], 0.0
    for name, k in (("cliques", 3), ("dimer", 3), ("p4", 10), ("singleton_chain", 6)):
        P = singleton_chain(k) if name == "singleton_chain" else nm.rig_setup(name, k)
        blocks = range(k) if k <= 6 else (1, 5, 9)
        for i in blocks:
            t0 = time.time()
            rd, _ = nm.n3x_upper_dense(P, i, a.M)
            t1 = time.time()
            rr, _ = nm.n3x_upper_reduced(P, i, a.M)
            t2 = time.time()
            for x in (P["s"], P["t"]):
                d, r = arb(nm.upper_q(rd[x])), arb(nm.upper_q(rr[x]))
                rel = abs(float((d - r) / d)) if nm.lower_q(d) > 0 else 0.0
                worst = max(worst, rel)
                rows.append(dict(family=name, k=k, block=i, x=x, M=a.M,
                                 log10_dense=f"{nm.lg10(d):.12f}", log10_reduced=f"{nm.lg10(r):.12f}",
                                 rel_diff=f"{rel:.3e}", sec_dense=f"{t1-t0:.2f}",
                                 sec_reduced=f"{t2-t1:.2f}"))
            print(f"{name} {k} block {i}: dense {t1-t0:.1f}s reduced {t2-t1:.1f}s")
    out = HERE.parent / "results" / "crosscheck.csv"
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"max relative difference dense vs reduced: {worst:.3e}  ->",
          "PASS" if worst < 1e-25 else "FAIL")
    sys.exit(0 if worst < 1e-25 else 1)
