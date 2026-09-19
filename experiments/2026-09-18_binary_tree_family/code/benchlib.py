"""Exact-diagonalization benchmarks next to the certificates of Tables 3-5 (R-B5).

The certificates are READ from the Arb run of the sibling experiment
(../2026-09-18_interval_certification/results/aggregate.csv); the diagonalization
here only measures how far they sit from the quantities they bound.

Usage: python run_benchmark.py [--smoke] [--dps 60] [--only fam:k,...] [--tag name]
"""
import argparse
import csv
import sys
import time
from fractions import Fraction

import mpmath as mp
import numpy as np

import bounds as BD
import logutil as LU
import spectral_mp as SP
from families import build, dense

CERT = LU.EXP.parent / "2026-09-18_interval_certification" / "results" / "aggregate.csv"
RES = LU.EXP / "results"
PROV = ("Addresses delegated review finding R-B5 (revision_response_trace.json): exact "
        "ITF S^2, time-averaged fidelity, F_max bracket and block-projector masses next "
        "to the Section 7 certificates, with certificate/true ratios. Calibration only: "
        "the certificates themselves never use these numbers.")


def mpq(x: Fraction):
    return mp.mpf(x.numerator) / x.denominator


def float64_S(G):
    H = np.array(dense(G, float))
    lam, Q = np.linalg.eigh(H)
    groups, cur = [], [0]
    for j in range(1, len(lam)):
        if lam[j] - lam[j - 1] < 1e-9:
            cur.append(j)
        else:
            groups.append(cur)
            cur = [j]
    groups.append(cur)
    s, t = G["s"], G["t"]
    return sum(abs(sum(Q[t, j] * Q[s, j] for j in g)) for g in groups)


def benchmark_instance(G, kappa, dps, cert_row=None, tol_exp=-40):
    mp.mp.dps = dps
    t0 = time.time()
    H = SP.to_mp(dense(G, mpq))
    lam, Q = SP.eig(H)
    res, orth = SP.residuals(H, lam, Q)
    groups, min_gap, within = SP.group(lam, mp.mpf(10) ** tol_exp)
    s, t = G["s"], G["t"]
    c = SP.coefficients(Q, groups, s, t)
    glam = [mp.fsum(lam[j] for j in g) / len(g) for g in groups]
    S, Fbar = SP.S_Fbar(c)
    Fs, fmeta = SP.F_sampled(c, glam)
    bd = BD.block_data(G)
    par = BD.derived(G, bd)
    cert = BD.certificates(G, par, mpq(Fraction(kappa)))
    intervals = [(lo - bd["delta"], hi + bd["delta"]) for lo, hi in bd["ends"]]
    ms, cnt = SP.block_masses(Q, lam, intervals, s)
    mt, _ = SP.block_masses(Q, lam, intervals, t)
    Sigma2 = mp.fsum((a * b) ** 2 for a, b in zip(cert["Ls"], cert["Lt"]))
    checks = dict(
        aglo_counts=all(cnt[i] == len(V) for i, V in enumerate(G["parts"])) and sum(cnt) == G["n"],
        S_le_Sigma1=S <= cert["Sigma1"],
        S_le_B=S ** 2 <= cert["B2"],
        Fbar_le_Sigma2=Fbar <= Sigma2,
        Fsamp_le_S2=Fs <= S ** 2 * (1 + mp.mpf(10) ** -10),
        masses_le_L=all(ms[i] <= cert["Ls"][i] and mt[i] <= cert["Lt"][i] for i in range(G["k"])),
    )
    out = dict(family=G["name"], k=G["k"], n=G["n"], D=cert["D"], kappa=str(kappa), dps=dps,
               S=S, ITF=S ** 2, Fbar=Fbar, F_samp=Fs,
               Sigma1=cert["Sigma1"], F_esp=cert["F_esp"], B2=cert["B2"], b_eta2=cert["b_eta2"],
               B_E2=cert["B_E2"], Sigma2=Sigma2, gamma_kappa=cert["gamma"],
               ratio_Sigma1_over_S=cert["Sigma1"] / S if S > 0 else mp.inf,
               ratio_B2_over_ITF=cert["B2"] / S ** 2 if S > 0 else mp.inf,
               ratio_Fesp_over_Fsamp=cert["F_esp"] / Fs if Fs > 0 else mp.inf,
               eig_residual=res, orth_error=orth, min_group_gap=min_gap, max_within_group=within,
               n_distinct=len(groups), fsamp_meta=fmeta, seconds=time.time() - t0)
    if cert_row is not None:
        up = mp.mpf(cert_row["Sigma1_upper_certified"])
        out["Sigma1_vs_arb_reldiff"] = abs(cert["Sigma1"] - up) / up
        checks["Sigma1_matches_arb"] = out["Sigma1_vs_arb_reldiff"] < mp.mpf(10) ** -9
    if S > 1e-10:
        out["S_float64"] = float64_S(G)
        out["S_float64_reldiff"] = abs(mp.mpf(out["S_float64"]) - S) / S
    masses = [dict(family=G["name"], k=G["k"], block=i, dist_s=BD.dists(G, s)[0][i],
                   dist_t=BD.dists(G, t)[0][i], mass_s=ms[i], L_s=cert["Ls"][i],
                   mass_t=mt[i], L_t=cert["Lt"][i],
                   ratio_s=cert["Ls"][i] / ms[i] if ms[i] > 0 else mp.inf,
                   ratio_t=cert["Lt"][i] / mt[i] if mt[i] > 0 else mp.inf)
              for i in range(G["k"])]
    return out, checks, masses


def adaptive_instance(G, kappa, dps0, cert_row=None, confirm=False):
    """Raise the working precision until S is resolved: S must exceed 10^-(dps-30).
    Then (if confirm or the precision had to grow) re-run at dps+20 and require the
    relative change in S to be < 1e-6.  Returns (out, checks, masses, dps_trail)."""
    dps, trail = dps0, []
    while True:
        out, checks, masses = benchmark_instance(G, kappa, dps, cert_row)
        trail.append((dps, mp.nstr(out["S"], 8), round(out["seconds"], 1)))
        if out["S"] > mp.mpf(10) ** (-(dps - 30)):
            break
        est = -int(mp.floor(mp.log10(out["S"]))) if out["S"] > 0 else dps
        dps = max(dps + 60, est + 60) if out["S"] > 0 else dps * 2
    if confirm or dps > dps0:
        S_ref = out["S"]
        out2, _, _ = benchmark_instance(G, kappa, dps + 20, cert_row)
        rel = abs(out2["S"] - S_ref) / S_ref
        trail.append((dps + 20, mp.nstr(out2["S"], 8), round(out2["seconds"], 1)))
        checks["precision_converged"] = rel < mp.mpf(10) ** -6
        out["S_precision_reldiff"] = rel
    out["dps_trail"] = str(trail)
    out["dps"] = dps
    return out, checks, masses


def fmt(v):
    if isinstance(v, (mp.mpf, float)):
        return mp.nstr(mp.mpf(v), 8, min_fixed=-3, max_fixed=4) if v != mp.inf else "inf"
    return v


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--dps", type=int, default=60)
    ap.add_argument("--only", default="")
    ap.add_argument("--tag", default="")
    a = ap.parse_args()
    LU.capture_env()
    hdr = LU.header("2026-09-18_exact_benchmarks", "numerical", ["numerical"], PROV,
                    "2026-09-18T00:00:00-03:00")
    rows = list(csv.DictReader(open(CERT)))
    if a.only:
        want = {tuple(x.split(":")) for x in a.only.split(",")}
        rows = [r for r in rows if (r["family"], r["k"]) in want]
    if a.smoke:
        rows = [r for r in rows if (r["family"], r["k"]) == ("cliques", "3")]
    RES.mkdir(exist_ok=True)
    agg, allm, bad = [], [], []
    for r in rows:
        ts, t0 = LU.now(), time.time()
        G = build(r["family"], int(r["k"]))
        out, checks, masses = adaptive_instance(
            G, Fraction(r["kappa_tab"]), a.dps, r,
            confirm=(r["family"], r["k"]) == ("p4", "40"))
        ok = all(checks.values())
        if not ok:
            bad.append((r["family"], r["k"], {k: v for k, v in checks.items() if not v}))
        agg.append({k: fmt(v) for k, v in out.items() if k != "fsamp_meta"} |
                   {"checks_all_pass": ok})
        allm += [{k: fmt(v) for k, v in m.items()} for m in masses]
        print(f"{r['family']:8s} k={r['k']:>3} n={G['n']:>4} S={mp.nstr(out['S'],6)} "
              f"Sigma1={mp.nstr(out['Sigma1'],6)} ratio={mp.nstr(out['ratio_Sigma1_over_S'],4)} "
              f"B2/ITF={mp.nstr(out['ratio_B2_over_ITF'],4)} checks={'ok' if ok else checks} "
              f"[{out['seconds']:.1f}s]", flush=True)
        LU.append_run(hdr, dict(
            timestamp_start=ts, timestamp_end=LU.now(), wallclock_seconds=round(time.time() - t0, 1),
            peak_memory_mb=round(LU.peak_mb(), 1), status="success" if ok else "failure", seed=None,
            parameters=dict(family=r["family"], k=int(r["k"]), kappa=r["kappa_tab"], dps=a.dps,
                            group_tol="1e-40", scale="smoke" if a.smoke else "full", tag=a.tag),
            inputs=dict(data_files=[str(CERT.relative_to(LU.ROOT))],
                        data_hashes={"interval_certification/aggregate.csv": LU.sha256_file(CERT)}),
            outputs=dict(scalar_results={k: (mp.nstr(v, 10) if isinstance(v, mp.mpf) else v)
                                         for k, v in out.items() if k not in ("fsamp_meta",)} |
                         {"F_samp_meta": out["fsamp_meta"], "checks": checks}),
            exit_code=0 if ok else 1, provenance_note=PROV,
            notes="seed=null: deterministic. F_samp is a lower estimate of F_max (sampled sup)."))
    name = "smoke" if a.smoke else (a.tag or "aggregate")
    for fname, data in ((f"{name}.csv" if name != "aggregate" else "aggregate.csv", agg),
                        (f"masses_{name}.csv" if name != "aggregate" else "masses.csv", allm)):
        fields = list(dict.fromkeys(k for row in data for k in row))   # union, first-seen order
        with open(RES / fname, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fields, restval="")
            w.writeheader()
            w.writerows(data)
    for b in bad:
        print("VIOLATION:", b)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
