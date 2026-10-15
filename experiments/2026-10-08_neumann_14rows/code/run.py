"""Run the Neumann-route certificates on the 14 rows of Table 3 (see ../plan.md).

  python run.py --smoke                 # K4 k=3, all bounds, N3x at M=200
  python run.py --tag final             # 14 rows, N3x at M=1000 and M=2000
  python run.py --tag prec256 --prec 256 --rows cliques:10,dimer:12,p4:40 --M 1000

Phase A (closed forms, one task per row): N1, N2, N2' (stated/path/ladder), N3, N3vN2.
Phase B (one task per row, block i and M): rigorous N3x upper sum (reduced system).
Every per-block bound is stored as the exact rational upper endpoint of its Arb ball
(results/raw_<tag>.json); aggregate.py-style sums are then exact rational arithmetic.
"""
import argparse
import json
import multiprocessing as mp
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
HERE = Path(__file__).resolve().parent
EXP = HERE.parent
sys.path.insert(0, str(HERE))

import flint                                    # noqa: E402
import logutil as lu                            # noqa: E402
import neumann as nm                            # noqa: E402

ROWS = [("cliques", k) for k in (3, 4, 5, 6, 8, 10)] + \
       [("dimer", k) for k in (3, 5, 7, 9, 12)] + [("p4", k) for k in (10, 20, 40)]
PROVENANCE = ("Phase 4 of opportunity #2 (LAA rejection, reasons 2 and 3): Neumann-route "
              "certificates N1, N2, N2', N3 and a rigorous N3x on the 14 rows of Table 3, "
              "against the exact S and block masses of 2026-09-18_exact_benchmarks. MR2 (rate "
              "share) and MR3 (N3x as a certificate). Calibrated by "
              "experiments/toy/2026-10-08_neumann_calibration.")
_P = {}


def _init(prec):
    flint.ctx.prec = prec


def _getP(name, k):
    if (name, k) not in _P:
        _P[(name, k)] = nm.rig_setup(name, k)
    return _P[(name, k)]


def _ser(d):
    return {f"{i},{x}": str(nm.upper_q(v)) for (i, x), v in d.items()}


def task(spec):
    t0, ts = time.time(), lu.now()
    try:
        if spec[0] == "closed":
            _, name, k = spec
            P = _getP(name, k)
            out, cert = nm.bounds_123(P)
            forms, info = nm.n2prime(P)
            out.update(forms)
            res = {c: _ser(d) for c, d in out.items()}
            extra = dict(cert=cert, info=info, g=str(P["g"]), delta=str(P["delta"]),
                         eta=str(P["eta"]), D=P["D"])
        else:
            _, name, k, i, M = spec
            P = _getP(name, k)
            r, nseg = nm.n3x_upper_reduced(P, i, M)
            res = {str(x): str(nm.upper_q(v)) for x, v in r.items()}
            extra = dict(segments=nseg)
        status = "success"
    except Exception as e:                       # logged, never retried (Rule 14)
        res, extra, status = {}, dict(error=repr(e)), "failure"
    return dict(spec=list(spec), res=res, extra=extra, status=status, t_start=ts,
                t_end=lu.now(), seconds=round(time.time() - t0, 3), peak_mb=lu.peak_mb())


def cost(spec):
    if spec[0] == "closed":
        return spec[2] ** 2
    return spec[2] ** 2 * spec[4] / 100


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--tag", default="final")
    ap.add_argument("--prec", type=int, default=128)
    ap.add_argument("--M", default="1000,2000")
    ap.add_argument("--rows", default="")
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--no-n3x", action="store_true")
    a = ap.parse_args()
    if a.smoke:
        a.tag, rows, Ms = "smoke", [("cliques", 3)], [200]
    else:
        rows = [(r.split(":")[0], int(r.split(":")[1])) for r in a.rows.split(",")] \
            if a.rows else ROWS
        Ms = [] if a.no_n3x else [int(m) for m in a.M.split(",")]
    flint.ctx.prec = a.prec
    info = lu.capture_env()
    hdr = lu.header("2026-10-08_neumann_14rows", "numerical", ["numerical", "symbolic"],
                    PROVENANCE, "2026-10-08T16:40:00-03:00")
    specs = [("closed", n, k) for n, k in rows]
    specs += [("n3x", n, k, i, M) for n, k in rows for M in Ms for i in range(k)]
    specs.sort(key=cost, reverse=True)
    raw_path = EXP / "results" / f"raw_{a.tag}.json"
    raw = json.loads(raw_path.read_text()) if raw_path.exists() else {}
    done = {tuple(r["spec"]) for r in raw.get("tasks", []) if r["status"] == "success"}
    todo = [s for s in specs if tuple(s) not in done]
    print(f"{len(specs)} tasks, {len(todo)} to run, {a.workers} workers, prec {a.prec}",
          flush=True)
    tasks = raw.get("tasks", [])
    t0 = time.time()
    with mp.get_context("fork").Pool(a.workers, initializer=_init, initargs=(a.prec,)) as pool:
        for n_done, r in enumerate(pool.imap_unordered(task, todo), 1):
            tasks = [t for t in tasks if t["spec"] != r["spec"]] + [r]
            raw = dict(tag=a.tag, prec=a.prec, Ms=Ms, rows=[list(x) for x in rows], tasks=tasks)
            raw_path.write_text(json.dumps(raw) + "\n")      # checkpoint after every task
            if r["status"] != "success" or n_done % 20 == 0 or r["spec"][0] == "closed":
                print(f"[{time.time()-t0:7.1f}s] {n_done}/{len(todo)} {r['spec']} "
                      f"{r['status']} {r['seconds']}s {r['extra'].get('error', '')}", flush=True)
    # one log entry per row (closed forms) and per (row, M) (N3x, aggregated over blocks)
    sha = {p: lu.sha256_file(p) for p in (
        HERE / "neumann.py", HERE / "run.py", HERE / "logutil.py",
        lu.ROOT / "experiments/toy/2026-10-08_neumann_calibration/code/run.py",
        lu.ROOT / "experiments/2026-09-18_exact_benchmarks/code/logutil.py",
        lu.ROOT / "experiments/2026-09-18_exact_benchmarks/code/families.py",
        lu.ROOT / "experiments/2026-09-18_interval_certification/code/certify.py")
        if p.exists()}  # public export: the calibration toy is not published
    inputs = dict(code_files=[str(p.relative_to(lu.ROOT)) for p in sha],
                  code_hashes={str(p.relative_to(lu.ROOT)): h for p, h in sha.items()})
    by = {}
    for r in tasks:
        s = r["spec"]
        key = ("closed", s[1], s[2]) if s[0] == "closed" else ("n3x", s[1], s[2], s[4])
        by.setdefault(key, []).append(r)
    for key, rs in sorted(by.items(), key=lambda kv: str(kv[0])):
        if not any(tuple(r["spec"]) in {tuple(t) for t in todo} for r in rs):
            continue                                   # logged by an earlier invocation
        ok = all(r["status"] == "success" for r in rs)
        params = dict(family=key[1], k=key[2], certificates=(
            "N1,N2,N2p_stated,N2p_path,N2p_ladder,N3,N3vN2" if key[0] == "closed"
            else f"N3x upper Riemann sum, M={key[3]}"), prec=a.prec, tag=a.tag,
            scale="smoke" if a.smoke else "canonical")
        if key[0] == "n3x":
            params["M"] = key[3]
            params["blocks"] = len(rs)
            params["segments_total"] = sum(r["extra"].get("segments", 0) for r in rs)
        lu.append_run(hdr, dict(
            timestamp_start=min(r["t_start"] for r in rs), timestamp_end=max(r["t_end"] for r in rs),
            wallclock_seconds=round(sum(r["seconds"] for r in rs), 2),
            peak_memory_mb=round(max(r["peak_mb"] for r in rs), 1),
            status="success" if ok else "failure", seed=None, parameters=params,
            inputs=inputs, outputs=dict(result_files=[str(raw_path.relative_to(EXP))]),
            exit_code=0 if ok else 1,
            notes="Deterministic interval arithmetic: no randomness, seed null. "
                  + ("wallclock_seconds is the sum of per-block task times (parallel pool)."
                     if key[0] == "n3x" else "")
                  + "; ".join(r["extra"].get("error", "") for r in rs if r["status"] != "success")))
    nfail = sum(r["status"] != "success" for r in tasks)
    print(f"done in {time.time()-t0:.1f}s; failures: {nfail}; raw -> {raw_path}")
    return 0 if nfail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
