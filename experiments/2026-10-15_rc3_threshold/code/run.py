"""Run the R-C3 threshold scan (see ../plan.md).

  python -I run.py --smoke                          # K4 k=3, published w, contour at M=200
  python -I run.py --only cliques:1/10:3,dimer:1/20:3,p4:3/10:10   # reproduction rows first
  python -I run.py                                  # all 121 rows (resumes from raw_final.json)

Tasks: ("closed", fam, w, k)      bridge, walk-sum, chain, Combes-Thomas (kappa grid 10^-4)
       ("contour", fam, w, k, i, M) contour certificate on block i
       ("exact", fam, w, k)       exact S and block masses (mpmath, adaptive dps)
Per-block bounds are stored as the exact rational upper endpoints of their Arb balls
(results/raw_<tag>.json); aggregate.py sums them in exact rational arithmetic.
"""
import argparse
import json
import multiprocessing as mp
import os
import sys
import time
from fractions import Fraction
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
HERE = Path(__file__).resolve().parent
EXP = HERE.parent
sys.path.insert(0, str(HERE))

import flint                                    # noqa: E402
import logutil as lu                            # noqa: E402  (this folder's, before rc3)
import rc3                                      # noqa: E402
import mpmath                                   # noqa: E402

assert Path(lu.__file__).resolve().parent == HERE, lu.__file__
PREC = 128
PROVENANCE = (
    "Closes the delegation of review card R-C3 (Phase 6 of opportunity #2, J. Phys. A; "
    "peer_review_reconciled.json and .sciflow/review_triage.json, accepted and delegated to "
    "experiment-runner on 2026-10-15). For each family of Tab:familias, the smallest chain "
    "length k from which each certificate (bridge, walk-sum, chain, contour, Combes-Thomas) "
    "guarantees Sigma_1 < sqrt2 - 1, hence F_Bloch < 2/3 by cor:fent, next to the exact S; "
    "on K4 and the dimer, also stronger bridges (q = 0.20-0.45) where the baseline fails and "
    "the contour certifies.")
_CACHE = {}


def _init(prec):
    flint.ctx.prec = prec


def _get(name, w, k):
    key = (name, w, k)
    if key not in _CACHE:
        _CACHE[key] = rc3.setup(name, k, Fraction(w))
    return _CACHE[key]


def _ser(d):
    return {f"{i},{x}": str(v) for (i, x), v in d.items()}


def task(spec):
    t0, ts = time.time(), lu.now()
    try:
        kind, name, w, k = spec[:4]
        G, P = _get(name, w, k)
        if kind == "closed":
            res, flags = rc3.closed(P, G)
            res = {c: _ser(d) for c, d in res.items()}
            extra = dict(flags=flags, g_lower=str(P["g"]), delta=str(P["delta"]),
                         eta_upper=str(P["eta"]), mu=str(G["mu"]), d=G["d"], n=G["n"])
        elif kind == "contour":
            i, M = spec[4], spec[5]
            d, nseg = rc3.contour(P, i, M)
            res, extra = _ser(d), dict(segments=nseg)
        else:
            ex = rc3.exact(G, P)
            res = dict(S=mpmath.nstr(ex["S"], 30), Fbar=mpmath.nstr(ex["Fbar"], 20),
                       mass_s=[mpmath.nstr(m, 25) for m in ex["ms"]],
                       mass_t=[mpmath.nstr(m, 25) for m in ex["mt"]],
                       certifies=ex["certifies"], margin_rel=ex["margin_rel"])
            extra = dict(dps=ex["dps"], dps_trail=ex["dps_trail"],
                         checks={c: bool(v) for c, v in ex["checks"].items()},
                         residual=mpmath.nstr(ex["res"], 3), orth=mpmath.nstr(ex["orth"], 3),
                         min_group_gap=mpmath.nstr(ex["min_gap"], 3),
                         n_distinct=ex["n_distinct"],
                         S_precision_reldiff=ex.get("S_precision_reldiff"))
        status = "success"
    except Exception as e:                       # logged, never retried
        res, extra, status = {}, dict(error=repr(e)), "failure"
    return dict(spec=list(spec), res=res, extra=extra, status=status, t_start=ts,
                t_end=lu.now(), seconds=round(time.time() - t0, 3), peak_mb=lu.peak_mb())


def cost(spec):
    k = spec[3]
    if spec[0] == "contour":
        return k * spec[5] / 10
    if spec[0] == "exact":
        return k ** 3
    return k * 50


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--tag", default="final")
    ap.add_argument("--M", type=int, default=rc3.M_CONTOUR)
    ap.add_argument("--only", default="", help="fam:w:k,... (w as p/q)")
    ap.add_argument("--workers", type=int, default=6)
    a = ap.parse_args()
    flint.ctx.prec = PREC
    rc3.selfcheck()
    allrows = [(f, str(w), k) for f, reg, w, q, k in rc3.rows()]
    if a.smoke:
        a.tag, a.M = "smoke", 200
        rows = [("cliques", "1/10", 3)]
    elif a.only:
        want = [tuple(x.split(":")) for x in a.only.split(",")]
        rows = [(f, str(Fraction(w)), int(k)) for f, w, k in want]
        assert all(r in allrows for r in rows), rows
    else:
        rows = allrows
    lu.capture_env()
    hdr = lu.header("2026-10-15_rc3_threshold", "numerical", ["numerical", "symbolic"],
                    PROVENANCE, "2026-10-15T21:00:00-03:00")
    specs = []
    for f, w, k in rows:
        specs += [("closed", f, w, k), ("exact", f, w, k)]
        specs += [("contour", f, w, k, i, a.M) for i in range(k)]
    specs.sort(key=cost, reverse=True)
    raw_path = EXP / "results" / f"raw_{a.tag}.json"
    raw = json.loads(raw_path.read_text()) if raw_path.exists() else {}
    done = {tuple(r["spec"]) for r in raw.get("tasks", []) if r["status"] == "success"}
    todo = [s for s in specs if tuple(s) not in done]
    print(f"{len(specs)} tasks, {len(todo)} to run, {a.workers} workers, prec {PREC}, "
          f"M {a.M}", flush=True)
    tasks = raw.get("tasks", [])
    t0 = time.time()
    with mp.get_context("fork").Pool(a.workers, initializer=_init, initargs=(PREC,)) as pool:
        for n_done, r in enumerate(pool.imap_unordered(task, todo), 1):
            tasks = [t for t in tasks if t["spec"] != r["spec"]] + [r]
            raw = dict(tag=a.tag, prec=PREC, M=a.M, tasks=tasks)
            raw_path.write_text(json.dumps(raw) + "\n")      # checkpoint after every task
            if r["status"] != "success" or r["spec"][0] != "contour" or n_done % 50 == 0:
                print(f"[{time.time()-t0:7.1f}s] {n_done}/{len(todo)} {r['spec']} "
                      f"{r['status']} {r['seconds']}s {r['extra'].get('error', '')}", flush=True)
    sha = {p: lu.sha256_file(p) for p in [HERE / "rc3.py", HERE / "run.py", HERE / "logutil.py"]
           + rc3.REUSED}
    inputs = dict(code_files=[str(p.relative_to(lu.ROOT)) for p in sha],
                  code_hashes={str(p.relative_to(lu.ROOT)): h for p, h in sha.items()})
    todo_set = {tuple(s) for s in todo}
    by = {}
    for r in tasks:
        s = r["spec"]
        by.setdefault((s[1], s[2], s[3]), []).append(r)
    for (f, w, k), rs in sorted(by.items(), key=lambda kv: (kv[0][0], Fraction(kv[0][1]), kv[0][2])):
        if not any(tuple(r["spec"]) in todo_set for r in rs):
            continue                                   # logged by an earlier invocation
        ok = all(r["status"] == "success" for r in rs) and len(rs) == k + 2
        reg = "published" if Fraction(w) == rc3.W_PUBLISHED[f] else "stronger"
        params = dict(family=f, regime=reg, w=w, k=k, M=a.M, prec_bits=PREC,
                      kappa_grid=f"j/{rc3.KAPPA_DEN}", dps0=rc3.DPS0, tag=a.tag,
                      certificates="bridge (N1), walk-sum (N2), chain (N2p_path), contour "
                                   "(n3x_upper_reduced), Combes-Thomas (certify.sigma1); exact S",
                      scale="smoke" if a.smoke else "canonical")
        ex = next((r for r in rs if r["spec"][0] == "exact" and r["status"] == "success"), None)
        lu.append_run(hdr, dict(
            timestamp_start=min(r["t_start"] for r in rs), timestamp_end=max(r["t_end"] for r in rs),
            wallclock_seconds=round(sum(r["seconds"] for r in rs), 2),
            peak_memory_mb=round(max(r["peak_mb"] for r in rs), 1),
            status="success" if ok else "failure", seed=None, parameters=params,
            inputs=inputs,
            outputs=dict(result_files=[str(raw_path.relative_to(EXP))],
                         scalar_results=dict(S=ex["res"]["S"], dps=ex["extra"]["dps"],
                                             exact_checks=ex["extra"]["checks"]) if ex else {}),
            exit_code=0 if ok else 1,
            notes="Deterministic ball and exact arithmetic: no randomness, seed null. "
                  "wallclock_seconds is the sum of the task times of this row (parallel pool). "
                  + "; ".join(r["extra"].get("error", "") for r in rs if r["status"] != "success")))
    nfail = sum(r["status"] != "success" for r in tasks)
    print(f"done in {time.time()-t0:.1f}s; failures: {nfail}; raw -> {raw_path}")
    return 0 if nfail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
