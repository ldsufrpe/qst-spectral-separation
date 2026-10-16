"""Robustness of the Combes-Thomas failures near the threshold (added after the run; see README).

For every CT row that fails with Sigma_1 within 5 % above sqrt2 - 1, evaluate the CT Sigma_1 in
float64 on the fine grid kappa = j/10^6, j = 1..999999, and report its minimum. A minimum above
the threshold shows that the failure is not an artefact of the grid of step 10^-4 used for the
certificate. (A float value below the threshold would call for a rigorous evaluation; it would
only help CT.) Also reports whether the CT factor is capped at eta on every block off the
communities of s and t, in which case Sigma_1 = 2 eta + (k - 2) eta^2 for every kappa.

  python -I ct_refine.py
"""
import csv
import math
import sys
import time
from fractions import Fraction
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
EXP = HERE.parent
sys.path.insert(0, str(HERE))
import logutil as lu                                  # noqa: E402  (this folder's, before rc3)
import rc3                                            # noqa: E402
import flint                                          # noqa: E402

THR = math.sqrt(2) - 1


def fl(x):
    return float(rc3.C.upper_q(x))


def main():
    flint.ctx.prec = 128
    t0, ts = time.time(), lu.now()
    agg = list(csv.DictReader(open(EXP / "results" / "aggregate.csv")))
    rows = [r for r in agg if r["certificate"] == "CT" and r["certifies"] == "False"
            and float(r["margin_rel"]) > -0.05]
    kap = np.arange(1, 10 ** 6) / 10 ** 6
    out = []
    for r in rows:
        f, w, k = r["family"], Fraction(r["w_exact"]), int(r["k"])
        G = rc3.build_w(f, k, w)
        par, geo = rc3.C.parameters(G), rc3.C.geometry(G)
        eta, ghat = fl(par["eta"]), fl(par["ghat"])
        gam = np.log1p(kap * ghat)
        c = 1 / (1 - kap)
        root = math.sqrt(1 - eta ** 2)
        tot = np.zeros_like(kap)
        capped = True
        for i, V in enumerate(G["parts"]):
            A = c * fl(par["rho"][i]) * math.sqrt(len(V)) / root
            Ls, Lt = [], []
            for x, d, LL in ((G["s"], geo["ds"][i], Ls), (G["t"], geo["dt"][i], Lt)):
                b = 1.0 if G["comm"][x] == i else eta
                v = np.minimum(b, A * np.exp(-gam * d))
                if G["comm"][x] != i and np.any(v < eta):
                    capped = False
                LL.append(v)
            tot += Ls[0] * Lt[0]
        j = int(np.argmin(tot))
        out.append(dict(family=f, w=r["w"], k=k, sigma_grid_1e4=r["sigma"],
                        sigma_min_fine_1e6=f"{tot[j]:.7f}", kappa_fine=f"{kap[j]:.6f}",
                        eta_capped_every_kappa=capped,
                        eta_floor=f"{2 * eta + (k - 2) * eta ** 2:.7f}",
                        still_fails=bool(tot[j] > THR)))
        print(out[-1], flush=True)
    path = EXP / "results" / "ct_refine.csv"
    with open(path, "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(out[0]))
        wr.writeheader()
        wr.writerows(out)
    ok = all(x["still_fails"] for x in out)
    hdr = lu.header("2026-10-15_rc3_threshold", "numerical", ["numerical", "symbolic"],
                    rc3_provenance(), "2026-10-15T21:00:00-03:00")
    lu.append_run(hdr, dict(
        timestamp_start=ts, timestamp_end=lu.now(), wallclock_seconds=round(time.time() - t0, 2),
        peak_memory_mb=round(lu.peak_mb(), 1), status="success", seed=None,
        parameters=dict(tag="ct_refine", check="CT failures within 5 % of sqrt2-1, kappa = j/10^6 "
                        "in float64", rows=len(out), scale="canonical"),
        inputs=dict(code_files=["experiments/2026-10-15_rc3_threshold/code/ct_refine.py"],
                    code_hashes={"ct_refine.py": lu.sha256_file(HERE / "ct_refine.py")},
                    data_files=["results/aggregate.csv"]),
        outputs=dict(result_files=["results/ct_refine.csv"],
                     scalar_results=dict(rows=len(out), all_still_fail=ok)),
        exit_code=0 if ok else 1,
        notes="Added after the run, not in plan.md: two CT failures at the onset lie 0.2 % above "
              "the threshold. Float check of the kappa grid only; the certificates are unchanged."))
    return 0 if ok else 1


def rc3_provenance():
    import json
    return json.loads((EXP / "experiment_log.json").read_text())["provenance_note"]


if __name__ == "__main__":
    sys.exit(main())
