"""Mode localization in weakly coupled substructures (R-B8): certificate vs exact modal
analysis for (i) band participation ||P_i e_s||^2 and (ii) time-uniform displacement
transmissibility sup_tau |[cos(tau sqrt K)]_{ts}|, both bounded through
S(s,t) = max_{||f||_inf <= 1} |f(K)_{ts}| <= Sigma1(s,t).

Usage: python run.py [--smoke] [--kmax 10]
"""
import argparse
import csv
import math
import sys
import time
from fractions import Fraction as Q

import mpmath as mp
import numpy as np
from scipy.optimize import minimize_scalar

import benchlib
import bounds as BD
import logutil as LU
import spectral_mp as SP
from families import dense
from structure import build_structure

RES = LU.EXP / "results"
PROV = ("Addresses delegated review finding R-B8 (revision_response_trace.json): a "
        "non-quantum structured-matrix example (stiffness matrix of weakly coupled "
        "substructures, vibration mode localization) where block-projector mass bounds "
        "are the natural quantity; certificate vs exact modal values.")


def cos_transmissibility(c, glam, max_samples=2_000_000):
    """Lower estimate of sup_{tau>=0} |sum_lambda c_lambda cos(tau sqrt(lambda))|."""
    S = mp.fsum(abs(x) for x in c)
    w = np.array([float(x / S) for x in c])
    om = np.sqrt(np.array([float(l) for l in glam]))
    spread = om.max() - om.min()
    dt = math.pi / (8 * spread)
    diffs = np.diff(np.sort(om))
    diffs = diffs[diffs > 1e-9]
    T = min(200 * 2 * math.pi / diffs.min(), max_samples * dt)
    best, top = 0.0, []
    taus_all = np.arange(0.0, T, dt)
    for a in range(0, len(taus_all), 200_000):
        taus = taus_all[a:a + 200_000]
        amp = np.abs(np.cos(np.outer(taus, om)) @ w)
        j = np.argpartition(amp, -5)[-5:]
        top += [(amp[i], taus[i]) for i in j]
    top.sort(reverse=True)
    f = lambda tau: -abs(np.cos(tau * om) @ w)
    for _, tau0 in top[:20]:
        r = minimize_scalar(f, bounds=(max(0.0, tau0 - dt), tau0 + dt), method="bounded",
                            options=dict(xatol=1e-12))
        best = max(best, -r.fun)
    return S * mp.mpf(best), dict(T=T, dt=dt)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--kmax", type=int, default=10)
    a = ap.parse_args()
    ks = [3] if a.smoke else list(range(3, a.kmax + 1))
    LU.capture_env()
    hdr = LU.header("2026-09-18_mode_localization", "numerical", ["numerical"], PROV,
                    "2026-09-18T00:00:00-03:00")
    RES.mkdir(exist_ok=True)
    rows, parts_rows, bad = [], [], 0
    for k in ks:
        t0, ts = time.time(), LU.now()
        G = build_structure(k)
        mp.mp.dps = 50
        bd = BD.block_data(G)
        par = BD.derived(G, bd)
        assert par["H2"], "H2 fails"
        kap = BD.best_kappa(G, par, den=1000)
        kq = Q(int(mp.nint(kap * 1000)), 1000)
        out, checks, masses = benchlib.adaptive_instance(G, kq, 60)
        # recompute spectral data at the accepted precision for the matrix-function checks
        mp.mp.dps = out["dps"]
        H = SP.to_mp(dense(G, lambda x: mp.mpf(x.numerator) / x.denominator))
        lam, Qm = SP.eig(H)
        assert lam[0] > 0, "K not positive definite"
        groups, _, _ = SP.group(lam, mp.mpf(10) ** -40)
        c = SP.coefficients(Qm, groups, G["s"], G["t"])
        glam = [mp.fsum(lam[j] for j in g) / len(g) for g in groups]
        S = mp.fsum(abs(x) for x in c)
        f_sign = mp.fsum((1 if x > 0 else -1) * x for x in c)       # f(lambda)=sign(c_lambda)
        checks["maxf_attained_by_sign"] = abs(f_sign - S) <= mp.mpf(10) ** (-(out["dps"] - 10)) * max(S, 1)
        cos_sup, cmeta = cos_transmissibility(c, glam)
        checks["cos_sup_le_S"] = cos_sup <= S * (1 + mp.mpf(10) ** -10)
        ok = all(checks.values())
        bad += not ok
        row = dict(k=k, n=G["n"], D=out["D"], g=benchlib.fmt(bd["g"]), W_max=benchlib.fmt(max(bd["W"])),
                   delta=benchlib.fmt(bd["delta"]), eta=benchlib.fmt(par["eta"]), kappa=str(kq),
                   omega2_min=benchlib.fmt(lam[0]), omega2_max=benchlib.fmt(lam[-1]),
                   S=benchlib.fmt(S), cos_transmissibility_sampled=benchlib.fmt(cos_sup),
                   Sigma1=benchlib.fmt(out["Sigma1"]), B=benchlib.fmt(mp.sqrt(out["B2"])),
                   b_eta=benchlib.fmt(mp.sqrt(out["b_eta2"])),
                   ratio_Sigma1_over_S=benchlib.fmt(out["ratio_Sigma1_over_S"]),
                   dps=out["dps"], checks_all_pass=ok)
        rows.append(row)
        for m in masses:
            parts_rows.append(dict(k=k, band=m["block"], dist_s=m["dist_s"],
                                   participation_s=benchlib.fmt(m["mass_s"] ** 2),
                                   bound_s=benchlib.fmt(min(1, m["L_s"] ** 2)),
                                   participation_t=benchlib.fmt(m["mass_t"] ** 2),
                                   bound_t=benchlib.fmt(min(1, m["L_t"] ** 2))))
        print(" ".join(f"{k_}={v}" for k_, v in row.items()), f"[{time.time()-t0:.1f}s]", flush=True)
        LU.append_run(hdr, dict(
            timestamp_start=ts, timestamp_end=LU.now(), wallclock_seconds=round(time.time() - t0, 1),
            peak_memory_mb=round(LU.peak_mb(), 1), status="success" if ok else "failure", seed=None,
            parameters=dict(k=k, kappa=str(kq), scale="smoke" if a.smoke else "full") | G["physical"],
            outputs=dict(scalar_results=row | {"checks": {k_: bool(v) for k_, v in checks.items()},
                                               "cos_sampling": cmeta}),
            exit_code=0 if ok else 1, provenance_note=PROV,
            notes="seed=null: deterministic. cos transmissibility is a sampled lower estimate."))
    tag = "smoke_" if a.smoke else ""
    for name, data in ((f"{tag}aggregate.csv", rows), (f"{tag}participation.csv", parts_rows)):
        with open(RES / name, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(data[0].keys()))
            w.writeheader()
            w.writerows(data)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
