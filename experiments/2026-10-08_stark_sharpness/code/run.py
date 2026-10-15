"""Stark-ladder sharpness test (see ../plan.md).

  python run.py --smoke              # Z at F/J=10, N<=12; one chain L=6, F/J=10
  python run.py --part Z             # F/J in {5,6,8.5,10,16,32,100}, N = 1..100
  python run.py --part chain         # 25 chains, N3x for L <= 21

Writes results/z_series.csv, results/chains.csv, results/chain_blocks.csv and appends to
experiment_log.json (one entry per F/J on Z, one per chain instance).
"""
import argparse
import csv
import math
import multiprocessing as mp
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
HERE = Path(__file__).resolve().parent
EXP = HERE.parent
sys.path.insert(0, str(HERE))

import flint                                            # noqa: E402
import mpmath as mpm                                    # noqa: E402
import numpy as np                                      # noqa: E402
from flint import arb, fmpq                             # noqa: E402

import logutil as lu                                    # noqa: E402
import stark as st                                      # noqa: E402
nm = st.nm

FZ = [fmpq(5), fmpq(6), fmpq(17, 2), fmpq(10), fmpq(16), fmpq(32), fmpq(100)]
FC = [fmpq(6), fmpq(17, 2), fmpq(10), fmpq(16), fmpq(32)]
LC = [6, 11, 21, 41, 61]
J = fmpq(1)
M_N3X = 2000
PROVENANCE = ("Phase 4 of opportunity #2, MR1: sharpness of the walk-sum certificate (N2, "
              "N2' forms) on the uniform Wannier-Stark ladder, against the exact Bessel solution "
              "on Z and exact high-precision diagonalization of finite chains. Calibrated by "
              "experiments/toy/2026-10-08_neumann_calibration (T3).")


def best_kappa_np(F, Jv, N, src_tgt, grid=999):
    """Grid-optimal kappa (step 1/(grid+1)) for the paper's Sigma1, in float (the bound is
    valid for every kappa; only the choice is made in float). src_tgt: list of
    (dist_to_s, dist_to_t) over blocks."""
    p = st.params(F, Jv)
    g, d, eta, Jf = float(p["g"]), float(p["delta"]), float(p["eta"]), float(p["J"])
    ghat = (g - 2 * d) / (4 * Jf)
    rho_i = 1 + 4 * d / (math.pi * (g - 2 * d))
    kap = np.arange(1, grid + 1) / (grid + 1)
    A = rho_i / (1 - kap) / math.sqrt(1 - eta ** 2)
    gam = np.log1p(kap * ghat)
    ds = np.array([a for a, _ in src_tgt], float)[:, None]
    dt = np.array([b for _, b in src_tgt], float)[:, None]
    Ls = np.where(ds == 0, 1.0, np.minimum(eta, A[None, :] * np.exp(-gam[None, :] * ds)))
    Lt = np.where(dt == 0, 1.0, np.minimum(eta, A[None, :] * np.exp(-gam[None, :] * dt)))
    tot = (Ls * Lt).sum(axis=0)
    return fmpq(int(np.argmin(tot)) + 1, grid + 1)


def lgs(x):
    return f"{nm.lg10(x):.6f}" if x is not None else ""


# ------------------------------------------------------------------ Z
def run_Z(F, Nmax, T=40):
    x = 2 * J / F
    t0, ts = time.time(), lu.now()
    rows, viol = [], 0
    for N in range(1, Nmax + 1):
        mN = st.z_exact_mass(N, x)
        S = st.z_exact_S(N, F, J, T)
        sup = arb(2 * x).bessel_j(N)
        kap = best_kappa_np(F, J, N, [(abs(i), abs(N - i)) for i in range(-30, N + 31)])
        rec = dict(F_over_J=str(F / J), N=N, log10_mass=lgs(mN), log10_S=lgs(S),
                   log10_supamp=lgs(sup), S_ge_supamp=nm.lower_q(S) >= nm.upper_q(sup),
                   kappa_sigma1=str(kap))
        for form in st.FORMS_Z:
            b = st.z_bound(form, N, F, J, kap)
            if b is None:
                continue
            Sg = st.z_sigma(form, N, F, J, kap, T)
            ok = nm.upper_q(b[1]) >= nm.lower_q(mN) and nm.upper_q(Sg) >= nm.lower_q(S)
            viol += not ok
            rec[f"{form}_log10_L"] = lgs(b[1])
            rec[f"{form}_log10_Sigma"] = lgs(Sg)
            rec[f"{form}_valid"] = ok
        rows.append(rec)
    # every block, both ends: bound >= exact mass for distances 0..Nmax+T
    for form in st.FORMS_Z:
        kap = fmpq(1, 2)
        for n in range(0, Nmax + T + 2):
            b = st.z_bound(form, n, F, J, kap)
            if b is not None and nm.upper_q(b[1]) < nm.lower_q(st.z_exact_mass(n, x)):
                viol += 1
    return rows, viol, ts, lu.now(), time.time() - t0


# ------------------------------------------------------------------ chains
def exact_adaptive(L, F, dps0):
    trail = []
    dps = dps0
    while True:
        ex = st.chain_exact(L, F, J, dps)
        trail.append((dps, mpm.nstr(ex["S"], 8)))
        if ex["S"] > mpm.mpf(10) ** (-(dps - 30)):
            break
        dps += 60
    ex2 = st.chain_exact(L, F, J, dps + 20)
    trail.append((dps + 20, mpm.nstr(ex2["S"], 8)))
    rel = abs(ex2["S"] - ex["S"]) / ex2["S"]
    return ex2, dps + 20, rel, trail


def chain_task(spec):
    L, Fs = spec
    F = fmpq(*[int(v) for v in Fs.split("/")]) if "/" in Fs else fmpq(int(Fs))
    flint.ctx.prec = 128
    t0, ts = time.time(), lu.now()
    try:
        P = st.chain_P(L, F, J)
        s, t = P["s"], P["t"]
        out, cert = nm.bounds_123(P)
        forms, info = nm.n2prime(P)
        out.update(forms)
        kap = best_kappa_np(F, J, L - 1, [(i, L - 1 - i) for i in range(L)])
        out["Sigma1"] = st.chain_sigma1(L, F, J, kap)
        n3x_ok = L <= 21
        if n3x_ok:
            out["N3x"] = {}
            for i in range(L):
                r, _ = nm.n3x_upper_reduced(P, i, M_N3X)
                for xx in (s, t):
                    out["N3x"][(i, xx)] = r[xx]
        SZ = st.z_exact_S(L - 1, F, J)
        dps0 = max(50, int(-nm.lg10(SZ)) + 40)
        ex, dps, rel, trail = exact_adaptive(L, F, dps0)
        mpm.mp.dps = dps
        S = arb(mpm.nstr(ex["S"], 30))
        mass = {kk: arb(mpm.nstr(v, 30)) for kk, v in ex["mass"].items()}
        tol = fmpq(999999, 1000000)
        rec = dict(L=L, N=L - 1, F_over_J=str(F / J), dps=dps, S_dps_reldiff=f"{float(rel):.2e}",
                   residual=mpm.nstr(ex["residual"], 3), orth=mpm.nstr(ex["orth"], 3),
                   log10_S=lgs(S), log10_S_Z=lgs(SZ), kappa_sigma1=str(kap),
                   log10_mass_end=lgs(mass[(L - 1, 0)]))
        blocks, viol, n2_le_path = [], 0, True
        for c, d in out.items():
            sg = sum((arb(nm.upper_q(d[(i, s)])) * arb(nm.upper_q(d[(i, t)])) for i in range(L)),
                     arb(0))
            v = sum(nm.upper_q(d[(i, xx)]) < nm.lower_q(mass[(i, xx)]) * tol
                    for i in range(L) for xx in (s, t))
            v += int(nm.upper_q(sg) < nm.lower_q(S) * tol)
            viol += v
            rec[f"{c}_log10_Sigma"] = lgs(sg)
            rec[f"{c}_log10_L_end"] = lgs(d[(L - 1, s)])
            rec[f"{c}_violations"] = v
            for i in range(L):
                blocks.append(dict(L=L, F_over_J=str(F / J), block=i, certificate=c,
                                   log10_L_s=lgs(d[(i, s)]), log10_L_t=lgs(d[(i, t)]),
                                   log10_mass_s=lgs(mass[(i, s)]), log10_mass_t=lgs(mass[(i, t)])))
        if "N2p_path" in out:
            n2_le_path = all(nm.upper_q(out["N2"][kk]) <= nm.upper_q(out["N2p_path"][kk])
                             for kk in out["N2"])
        rec["N2_le_N2p_path"] = n2_le_path
        rec["rho_certified"] = cert["N2"] and cert["N3"]
        rec["violations_total"] = viol
        status = "success" if viol == 0 and rel < 1e-6 else "failure"
        return dict(spec=list(spec), rec=rec, blocks=blocks, status=status, trail=str(trail),
                    t_start=ts, t_end=lu.now(), seconds=round(time.time() - t0, 2),
                    peak_mb=lu.peak_mb(), error="")
    except Exception as e:
        return dict(spec=list(spec), rec={}, blocks=[], status="failure", trail="", t_start=ts,
                    t_end=lu.now(), seconds=round(time.time() - t0, 2), peak_mb=lu.peak_mb(),
                    error=repr(e))


def write(path, rows):
    keys = []
    for r in rows:
        keys += [k for k in r if k not in keys]
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--part", default="Z")
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--workers", type=int, default=6)
    a = ap.parse_args()
    flint.ctx.prec = 128
    lu.capture_env()
    hdr = lu.header("2026-10-08_stark_sharpness", "numerical", ["numerical", "symbolic"],
                    PROVENANCE, "2026-10-08T16:45:00-03:00")
    sha = {str(p.relative_to(lu.ROOT)): lu.sha256_file(p) for p in (
        HERE / "stark.py", HERE / "run.py", HERE / "logutil.py",
        lu.ROOT / "experiments/2026-10-08_neumann_14rows/code/neumann.py",
        lu.ROOT / "experiments/toy/2026-10-08_neumann_calibration/code/run.py",
        lu.ROOT / "experiments/2026-09-18_exact_benchmarks/code/spectral_mp.py")
        if p.exists()}  # public export: the calibration toy is not published
    inputs = dict(code_files=list(sha), code_hashes=sha)
    tag = "smoke_" if a.smoke else ""
    nfail = 0
    if a.part == "Z" or a.smoke:
        allrows = []
        for F in ([fmpq(10)] if a.smoke else FZ):
            rows, viol, ts, te, sec = run_Z(F, 12 if a.smoke else 100)
            allrows += rows
            ok = viol == 0 and all(r["S_ge_supamp"] for r in rows)
            nfail += not ok
            lu.append_run(hdr, dict(
                timestamp_start=ts, timestamp_end=te, wallclock_seconds=round(sec, 2),
                peak_memory_mb=round(lu.peak_mb(), 1), status="success" if ok else "failure",
                seed=None, parameters=dict(part="Z", F_over_J=str(F / J), J=1,
                                           N_max=len(rows), forms=",".join(st.FORMS_Z),
                                           prec=128, scale="smoke" if a.smoke else "canonical"),
                inputs=inputs, outputs=dict(result_files=[f"results/{tag}z_series.csv"],
                                            scalar_results=dict(violations=viol)),
                exit_code=0 if ok else 1,
                notes="Deterministic: seed null. Sigma1 kappa grid-optimal per N (float choice)."))
            print(f"Z F/J={F}: {len(rows)} N, violations {viol}, {sec:.1f}s", flush=True)
        write(EXP / "results" / f"{tag}z_series.csv", allrows)
    if a.part == "chain" or a.smoke:
        specs = [(6, "10")] if a.smoke else [(L, str(F)) for L in LC for F in FC]
        specs.sort(key=lambda s: -s[0])
        recs, blocks = [], []
        with mp.get_context("fork").Pool(a.workers) as pool:
            for r in pool.imap_unordered(chain_task, specs):
                nfail += r["status"] != "success"
                if r["rec"]:
                    recs.append(r["rec"])
                    blocks += r["blocks"]
                print(f"chain L={r['spec'][0]} F/J={r['spec'][1]}: {r['status']} {r['seconds']}s "
                      f"{r['error']}", flush=True)
                lu.append_run(hdr, dict(
                    timestamp_start=r["t_start"], timestamp_end=r["t_end"],
                    wallclock_seconds=r["seconds"], peak_memory_mb=round(r["peak_mb"], 1),
                    status=r["status"], seed=None,
                    parameters=dict(part="chain", L=r["spec"][0], F_over_J=r["spec"][1], J=1,
                                    n3x_M=M_N3X if r["spec"][0] <= 21 else None, prec=128,
                                    scale="smoke" if a.smoke else "canonical"),
                    inputs=inputs, outputs=dict(
                        result_files=[f"results/{tag}chains.csv", f"results/{tag}chain_blocks.csv"],
                        scalar_results=dict(violations=r["rec"].get("violations_total"),
                                            dps=r["rec"].get("dps"))),
                    exit_code=0 if r["status"] == "success" else 1,
                    notes=("Deterministic: seed null. dps trail " + r["trail"] + ". " + r["error"])))
        recs.sort(key=lambda r: (r["L"], float(fmpq(*[int(v) for v in r["F_over_J"].split("/")])
                                               if "/" in r["F_over_J"] else int(r["F_over_J"]))))
        write(EXP / "results" / f"{tag}chains.csv", recs)
        write(EXP / "results" / f"{tag}chain_blocks.csv", blocks)
    print(f"failures: {nfail}")
    return 0 if nfail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
