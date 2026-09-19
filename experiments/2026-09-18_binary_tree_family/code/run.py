"""Binary tree of K4 modules: Theorem C through (GR), beyond the chain corollary (R-B4).

Usage: python run.py [--rmax 12] [--rexact 5] [--smoke]
  Phase A  (all r <= rmax, no diagonalization): hypotheses, exact C(a) and C_theta.
  Phase B  choice of (kappa, beta, a) minimizing the distance D* at which the
           Theorem C bound  C e^{-2 beta D}  drops below 1e-2.
  Phase C  per-r certificates: Sigma1 (sharp), C1^unif e^{-gamma D} (eq:C1unif),
           Theorem C  min{1, C e^{-2 beta D}}.
  Phase D  (r <= rexact) exact S by multiprecision diagonalization, for calibration.
"""
import argparse
import csv
import math
import sys
import time
from fractions import Fraction as Q

import mpmath as mp
import numpy as np

import benchlib
import logutil as LU
from families import bfs, bridges_form_matching
from radial import radial_constants
from tree import build_tree

RES = LU.EXP / "results"
PROV = ("Addresses delegated review finding R-B4 (revision_response_trace.json): a "
        "non-chain family (binary tree of K4 modules) satisfying (GR) with a > 0, exercising "
        "Theorem C beyond Corollary cadeiasradiais, with explicit radial constants C, a.")

W, G_GAP, DELTA, D_MU, D0, M0 = 4.0, 20.0, 0.1, 4.0, 1, 4
ETA = DELTA / (G_GAP - DELTA)
GHAT = (G_GAP - 2 * DELTA) / (2 * D_MU)
RHO = 1 + 2 * (W + 2 * DELTA) / (math.pi * (G_GAP - 2 * DELTA))
A_GRID = [round(0.30 + 0.01 * j, 2) for j in range(51)] + [math.log(2) / 2]
TH_GRID = [round(0.30 + 0.01 * j, 2) for j in range(131)]
KAPPAS = [round(0.50 + 0.05 * j, 2) for j in range(10)]


def consts(kappa):
    gam = math.log(1 + kappa * GHAT)
    c = 1 / (1 - kappa)
    A0 = c * RHO * math.sqrt(M0) / math.sqrt(1 - ETA ** 2)
    return gam, c, A0


def check_hypotheses(G):
    assert bridges_form_matching(G), "bridges not a matching"
    assert G["d"] == 4 and G["mu"] == 1
    diam = 0
    for V in G["parts"][: min(len(G["parts"]), 64)]:
        for u in V:
            du = bfs(G["adj"], u) if G["n"] <= 4096 else None
            if du is not None:
                diam = max(diam, max(du[v] for v in V))
    D = bfs(G["adj"], G["s"])[G["t"]]
    assert D == 2 * G["r"] + 1, (D, G["r"])
    return dict(D=D, D0_sampled=diam)


def sigma1_float(G, kappa, ds, dt):
    gam, c, _ = consts(kappa)
    A = c * RHO * math.sqrt(M0) / math.sqrt(1 - ETA ** 2)
    si, ti = G["comm"][G["s"]], G["comm"][G["t"]]
    ds, dt = np.asarray(ds, float), np.asarray(dt, float)
    Ls = np.minimum(np.where(np.arange(G["k"]) == si, 1.0, ETA), A * np.exp(-gam * ds))
    Lt = np.minimum(np.where(np.arange(G["k"]) == ti, 1.0, ETA), A * np.exp(-gam * dt))
    return float(np.sum(Ls * Lt))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rmax", type=int, default=12)
    ap.add_argument("--rexact", type=int, default=5)
    ap.add_argument("--smoke", action="store_true")
    a = ap.parse_args()
    if a.smoke:
        a.rmax, a.rexact = 3, 2
    LU.capture_env()
    hdr = LU.header("2026-09-18_binary_tree_family", "numerical", ["numerical", "networks"],
                    PROV, "2026-09-18T00:00:00-03:00")
    RES.mkdir(exist_ok=True)
    tag = "smoke_" if a.smoke else ""

    # ---------------- Phase A
    per_r, radial_rows = {}, []
    for r in range(1, a.rmax + 1):
        t0, ts = time.time(), LU.now()
        G = build_tree(r)
        hyp = check_hypotheses(G)
        Ca, Ct = radial_constants(G, A_GRID, TH_GRID)
        per_r[r] = dict(G=G, hyp=hyp, Ca=Ca, Ct=Ct)
        for aa in A_GRID:
            radial_rows.append(dict(r=r, k=G["k"], quantity="C(a)", param=round(aa, 6),
                                    value=Ca[aa][0], argmax_vertex=Ca[aa][1]))
        for th in TH_GRID:
            radial_rows.append(dict(r=r, k=G["k"], quantity="C_theta", param=th,
                                    value=Ct[th][0], argmax_vertex=Ct[th][1]))
        print(f"[A] r={r:2d} k={G['k']:5d} n={G['n']:6d} D={hyp['D']} "
              f"C(ln2/2)={Ca[A_GRID[-1]][0]:.4f} C(0.40)={Ca[0.4][0]:.4f} "
              f"C_th(0.8)={Ct[0.8][0]:.4f} [{time.time()-t0:.1f}s]", flush=True)
        LU.append_run(hdr, dict(
            timestamp_start=ts, timestamp_end=LU.now(), wallclock_seconds=round(time.time() - t0, 1),
            peak_memory_mb=round(LU.peak_mb(), 1), status="success", seed=None,
            parameters=dict(phase="A_radial", r=r, k=G["k"], n=G["n"], scale="smoke" if a.smoke else "full"),
            outputs=dict(scalar_results=dict(D=hyp["D"], C_a_ln2_over_2=Ca[A_GRID[-1]][0],
                                             C_a={str(k_): v[0] for k_, v in Ca.items()},
                                             C_theta={str(k_): v[0] for k_, v in Ct.items()})),
            exit_code=0, provenance_note=PROV, notes="seed=null: deterministic BFS counts."))
    with open(RES / f"{tag}radial_constants.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(radial_rows[0].keys()))
        w.writeheader()
        w.writerows(radial_rows)

    # (GR) constant used in Theorem C: ANALYTIC, uniform in r (see README):
    #   dist(x,V_i) >= 2 d_T(j,i) - 1 and #{modules at tree distance rho} <= 3*2^(rho-1)
    #   => N_x(R) <= 3*2^((R+1)/2) - 2 <= 3*sqrt(2) * 2^(R/2)  => (GR) with a >= ln2/2, C = 3 sqrt 2.
    # The exact BFS values C(a) must never exceed it (checked for every r and a >= ln2/2).
    A_MIN = math.log(2) / 2
    C_AN = 3 * math.sqrt(2)
    viol = [(r, aa, per_r[r]["Ca"][aa][0]) for r in per_r for aa in A_GRID
            if aa >= A_MIN - 1e-12 and per_r[r]["Ca"][aa][0] > C_AN * (1 + 1e-12)]
    assert not viol, f"computed C(a) exceeds analytic 3*sqrt2: {viol[:3]}"
    Ctfam = {th: max(per_r[r]["Ct"][th][0] for r in per_r) for th in TH_GRID}
    a_grid_B = sorted(aa for aa in A_GRID if aa >= A_MIN - 1e-12)

    # ---------------- Phase B
    best = None
    for kap in KAPPAS:
        gam, c, A0 = consts(kap)
        for aa in a_grid_B:
            for jb in range(1, 400):
                beta = 0.005 * jb
                alpha = gam - beta
                if alpha <= aa + 1e-9:
                    break
                Cc = A0 ** 4 * math.exp(2 * beta * D0) * (C_AN * alpha / (alpha - aa)) ** 2
                Dstar = math.log(100 * Cc) / (2 * beta)
                if best is None or Dstar < best["Dstar"]:
                    best = dict(kappa=kap, a=aa, beta=beta, alpha=alpha, gamma=gam, A0=A0,
                                C_GR=C_AN, calC=Cc, Dstar=Dstar)
    best["C_GR_computed_max_over_r"] = max(per_r[r]["Ca"][best["a"]][0] for r in per_r)
    # exact-C_alpha variant of the same constant (tighter, same theorem: eq:radialF)
    th = min(TH_GRID, key=lambda x: abs(x - best["alpha"]) if x <= best["alpha"] else 9)
    best["C_alpha_exact_family"] = Ctfam[th]
    best["theta_used_for_C_alpha"] = th
    best["calC_exact_Calpha"] = best["A0"] ** 4 * math.exp(2 * best["beta"] * D0) * Ctfam[th] ** 2
    print("[B] choice:", {k_: (round(v, 6) if isinstance(v, float) else v) for k_, v in best.items()},
          flush=True)

    # ---------------- Phase C + D
    kap, beta, gam = best["kappa"], best["beta"], best["gamma"]
    _, c, A0 = consts(kap)
    rows, bad = [], 0
    for r, dat in per_r.items():
        t0, ts = time.time(), LU.now()
        G = dat["G"]
        ds_all = bfs(G["adj"], G["s"])
        dt_all = bfs(G["adj"], G["t"])
        ds = [min(ds_all[v] for v in V) for V in G["parts"]]
        dt = [min(dt_all[v] for v in V) for V in G["parts"]]
        D = dat["hyp"]["D"]
        kbest = min((j / 1000 for j in range(1, 1000)), key=lambda x: sigma1_float(G, x, ds, dt))
        S1 = sigma1_float(G, kbest, ds, dt)
        g2, c2, _ = consts(kbest)
        C1u = ((2 * c2 * RHO * math.sqrt(M0) / math.sqrt(1 - ETA ** 2)
                + c2 ** 2 * RHO ** 2 * G["k"] * M0 / (1 - ETA ** 2)) * math.exp(g2 * D0))
        thmC = min(1.0, best["calC"] * math.exp(-2 * beta * D))
        thmC_exact = min(1.0, best["calC_exact_Calpha"] * math.exp(-2 * beta * D))
        row = dict(r=r, k=G["k"], n=G["n"], D=D, kappa_sigma1=kbest,
                   Sigma1_sq=min(1.0, S1 ** 2), C1unif=C1u,
                   C1unif_bound_sq=min(1.0, (C1u * math.exp(-g2 * D)) ** 2),
                   thmC_bound=thmC, thmC_bound_exact_Calpha=thmC_exact,
                   C_GR_at_a=dat["Ca"][best["a"]][0], ITF="", ratio_Sigma1sq_over_ITF="",
                   ratio_thmC_over_ITF="", dps="")
        checks = {}
        if r <= a.rexact:
            out, chk, _ = benchlib.adaptive_instance(G, Q(round(kbest * 1000), 1000), 60)
            checks = dict(chk)
            ITF = out["ITF"]
            checks["ITF_le_thmC"] = ITF <= thmC
            checks["ITF_le_C1unif"] = ITF <= row["C1unif_bound_sq"]
            row.update(ITF=mp.nstr(ITF, 8), ratio_Sigma1sq_over_ITF=mp.nstr(row["Sigma1_sq"] / ITF, 6),
                       ratio_thmC_over_ITF=mp.nstr(thmC / ITF, 6), dps=out["dps"])
        ok = all(checks.values()) if checks else True
        bad += not ok
        rows.append(row | {"checks_all_pass": ok})
        print(f"[C] r={r:2d} D={D:2d} Sigma1^2={row['Sigma1_sq']:.3e} C1unif^2={row['C1unif_bound_sq']:.3e} "
              f"thmC={thmC:.3e} ITF={row['ITF'] or '-'} ok={ok} [{time.time()-t0:.1f}s]", flush=True)
        LU.append_run(hdr, dict(
            timestamp_start=ts, timestamp_end=LU.now(), wallclock_seconds=round(time.time() - t0, 1),
            peak_memory_mb=round(LU.peak_mb(), 1), status="success" if ok else "failure", seed=None,
            parameters=dict(phase="C_D_certificates", r=r, kappa_thmC=kap, beta=beta, a=best["a"],
                            kappa_sigma1=kbest, exact=r <= a.rexact),
            outputs=dict(scalar_results={k_: (str(v) if not isinstance(v, (int, float, str)) else v)
                                         for k_, v in row.items()} | {"checks": checks}),
            exit_code=0 if ok else 1, provenance_note=PROV,
            notes="seed=null: deterministic."))
    with open(RES / f"{tag}aggregate.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    with open(RES / f"{tag}theoremC_choice.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(best.keys()) + ["eta", "ghat", "rho", "D0", "M0"])
        w.writeheader()
        w.writerow(best | dict(eta=ETA, ghat=GHAT, rho=RHO, D0=D0, M0=M0))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
