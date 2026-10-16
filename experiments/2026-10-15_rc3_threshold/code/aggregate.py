"""Aggregate raw_<tag>.json into results/{aggregate,thresholds,reproduction,probe_agreement}.csv
and results/viz_schema.json, with the checks of plan.md (verification steps 3-6); with --final,
write the summary into experiment_log.json.

  python -I aggregate.py [--tag final] [--final]

Per-block bounds are exact rational upper endpoints, so every Sigma_1 is exact rational arithmetic
on certified upper bounds, and its verdict (Sigma_1 + 1)^2 < 2 is exact. S and the block masses
are mpmath values (calibration only); their verdict is numerical, with its margin.
"""
import argparse
import csv
import json
import math
import re
import sys
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
EXP = HERE.parent
sys.path.insert(0, str(HERE))
import logutil as lu                                  # noqa: E402  (this folder's, before rc3)
import rc3                                            # noqa: E402
import flint                                          # noqa: E402
import mpmath as mp                                   # noqa: E402
from flint import fmpq                                # noqa: E402

ROOT = rc3.ROOT
ROWS14 = ROOT / "experiments/2026-10-08_neumann_14rows/results/aggregate.csv"
CTREF = ROOT / "experiments/2026-09-18_interval_certification/results/aggregate.csv"
# public export: final.csv is byte-identical to the published aggregate.csv
SREF = ROOT / "experiments/2026-09-18_exact_benchmarks/results/aggregate.csv"
PROBE = EXP / "exploration/preplan_probe/probe_output.txt"
REF14 = {"bridge": "N1", "walk": "N2", "chain": "N2p_path", "contour": "N3x_M2000", "CT": "Sigma1"}
REPRO_ROWS = (("cliques", 3), ("dimer", 3), ("p4", 10))
ALL = ("exact",) + rc3.CERTS
THR = math.sqrt(2) - 1


def rd(path):
    with open(path) as f:
        return list(csv.DictReader(f))


def write(path, rows):
    keys = []
    for r in rows:
        for kk in r:
            if kk not in keys:
                keys.append(kk)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys, restval="")
        w.writeheader()
        w.writerows(rows)


def fq(s):
    p = s.split("/")
    return fmpq(int(p[0]), int(p[1])) if len(p) == 2 else fmpq(int(p[0]))


def blk(d):
    return {tuple(int(v) for v in kk.split(",")): fq(val) for kk, val in d.items()}


def wdec(w):
    return f"{float(Fraction(w)):g}"


def load(tag):
    raw = json.loads((EXP / "results" / f"raw_{tag}.json").read_text())
    R = {}
    for t in raw["tasks"]:
        if t["status"] != "success":
            continue
        s = t["spec"]
        key = (s[1], s[2], s[3])
        r = R.setdefault(key, dict(cert={}, contour={}, contour_blocks=set()))
        if s[0] == "closed":
            r["cert"].update({c: blk(d) for c, d in t["res"].items()})
            r["closed_extra"] = t["extra"]
        elif s[0] == "contour":
            r["contour"].update(blk(t["res"]))
            r["contour_blocks"].add(s[4])
            r["M"] = s[5]
        else:
            r["exact"], r["exact_extra"] = t["res"], t["extra"]
    complete = {}
    for (f, w, k), r in R.items():
        if "closed_extra" in r and "exact" in r and r["contour_blocks"] == set(range(k)):
            r["cert"]["contour"] = r["contour"]
            complete[(f, w, k)] = r
    return raw, complete


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="final")
    ap.add_argument("--final", action="store_true", help="all 121 rows required; write summary")
    a = ap.parse_args()
    flint.ctx.prec = rc3.PREC if hasattr(rc3, "PREC") else 128
    mp.mp.dps = 50
    raw, R = load(a.tag)
    meta = {(f, str(w), k): (reg, q) for f, reg, w, q, k in rc3.rows()}
    problems = []
    if a.final and set(R) != set(meta):
        problems.append(f"rows complete {len(R)}/{len(meta)}")
    agg, verdicts = [], {}
    for (f, w, k), r in sorted(R.items(), key=lambda kv: (kv[0][0], Fraction(kv[0][1]), kv[0][2])):
        reg, qt = meta[(f, w, k)]
        ce, ee = r["closed_extra"], r["exact_extra"]
        G, P = rc3.setup(f, k, Fraction(w))            # cheap; needed for s, t and k
        s, t = P["s"], P["t"]
        S = mp.mpf(r["exact"]["S"])
        ms = [mp.mpf(x) for x in r["exact"]["mass_s"]]
        mt = [mp.mpf(x) for x in r["exact"]["mass_t"]]
        if not all(ee["checks"].values()):
            problems.append(f"{f} w={w} k={k}: exact checks {ee['checks']}")
        if not ce["flags"]["rho_walk_certified"] or not ce["flags"]["chain_form_ok"]:
            problems.append(f"{f} w={w} k={k}: flags {ce['flags']}")
        base = dict(family=f, regime=reg, w=wdec(w), w_exact=w,
                    q_target=f"{float(qt):.2f}" if qt is not None else "",
                    q=f"{2 * float(Fraction(w)) / rc3.G_CLOSED[f]:.4f}",
                    eta=f"{float(fq(ce['eta_upper'])):.4f}", mu=ce["mu"], k=k, n=ce["n"])
        rec = dict(base, certificate="exact", log10_sigma=f"{float(mp.log10(S)):.6f}",
                   sigma=mp.nstr(S, 7), certifies=r["exact"]["certifies"],
                   margin_rel=f"{float((THR - S) / THR):.6f}", violations=0, dps=ee["dps"])
        agg.append(rec)
        verdicts.setdefault((f, w, "exact"), {})[k] = r["exact"]["certifies"]
        for c in rc3.CERTS:
            d = r["cert"].get(c)
            if d is None:
                problems.append(f"{f} w={w} k={k}: {c} missing")
                continue
            sig = rc3.sigma(P, d)
            viol = 0
            for i in range(k):
                for x, m in ((s, ms[i]), (t, mt[i])):
                    if mp.mpf(int(d[(i, x)].p)) / int(d[(i, x)].q) < m * (1 - mp.mpf(10) ** -6):
                        viol += 1
            sigf = mp.mpf(int(sig.p)) / int(sig.q)
            if sigf < S * (1 - mp.mpf(10) ** -6):
                viol += 1
            ok = rc3.certifies(sig)
            rec = dict(base, certificate=c, log10_sigma=f"{rc3.log10q(sig):.6f}",
                       sigma=mp.nstr(sigf, 7), certifies=ok,
                       margin_rel=f"{float((THR - sigf) / THR):.6f}", violations=viol)
            if c == "CT":
                rec["kappa"] = ce["flags"]["CT"]["kappa"]
                rec["kappa_argmin_size"] = ce["flags"]["CT"]["kappa_argmin_size"]
            if c == "contour":
                rec["M"] = r["M"]
            agg.append(rec)
            verdicts.setdefault((f, w, c), {})[k] = ok
            if viol:
                problems.append(f"{f} w={w} k={k} {c}: {viol} violations")
    write(EXP / "results" / "aggregate.csv", agg)

    # ---- reproduction (verification step 3)
    rows14 = {(x["family"], int(x["k"]), x["certificate"]): x for x in rd(ROWS14)}
    ctref = {(x["family"], int(x["k"])): x for x in rd(CTREF)}
    sref = {(x["family"], int(x["k"])): x for x in rd(SREF)}
    rep = []
    for f, k in REPRO_ROWS:
        w = str(rc3.W_PUBLISHED[f])
        mine = {x["certificate"]: x for x in agg if (x["family"], x["w_exact"], x["k"]) == (f, w, k)}
        if not mine:
            continue
        for c in rc3.CERTS:
            ours = float(mine[c]["log10_sigma"])
            if c == "CT":
                ref = math.log10(float(ctref[(f, k)]["Sigma1_upper_certified"]))
                src, tol = "interval_certification Sigma1_upper_certified (12 s.f.)", 1e-6 / math.log(10)
            else:
                ref = float(rows14[(f, k, REF14[c])]["log10_sigma"])
                src, tol = f"neumann_14rows {REF14[c]} log10_sigma (6 decimals)", 2e-6
            rep.append(dict(family=f, k=k, certificate=c, ours_log10=f"{ours:.7f}",
                            ref_log10=f"{ref:.7f}", diff=f"{ours - ref:+.2e}", tol=f"{tol:.1e}",
                            passes=abs(ours - ref) <= tol, reference=src))
        S_ours = mp.mpf(mine["exact"]["sigma"])
        S_ref = mp.mpf(sref[(f, k)]["S"])
        rel = float(abs(S_ours - S_ref) / S_ref)
        rep.append(dict(family=f, k=k, certificate="exact", ours_log10=mine["exact"]["log10_sigma"],
                        ref_log10=f"{float(mp.log10(S_ref)):.7f}", diff=f"{rel:+.2e} rel",
                        tol="1.0e-06 rel", passes=rel <= 1e-6,
                        reference="exact_benchmarks final.csv S (8 s.f.)"))
    write(EXP / "results" / "reproduction.csv", rep)
    for x in rep:
        if not x["passes"]:
            problems.append(f"reproduction {x['family']} k={x['k']} {x['certificate']}: {x['diff']}")

    # ---- probe agreement (verification step 5)
    probe, cur = {}, None
    for line in PROBE.read_text().splitlines():
        m = re.match(r"--- (\w+) w=([\d.]+)", line)
        if m:
            cur = (m.group(1), Fraction(m.group(2)))
            continue
        m = re.match(r"k=\s*(\d+) eta=\S+ (.*)", line)
        if m and cur:
            for name, val in re.findall(r"(\w+)=([\d.e+-]+)", m.group(2)):
                probe[(cur[0], cur[1], int(m.group(1)), name)] = float(val)
    pmap = {"exact": "S", "contour": "contour", "chain": "chain", "walk": "walk",
            "bridge": "bridge", "CT": "CT"}
    pa = []
    for x in agg:
        key = (x["family"], Fraction(x["w_exact"]), x["k"], pmap[x["certificate"]])
        if key not in probe:
            continue
        pv = probe[key]
        sv = float(x["sigma"])
        clear = abs(pv - THR) / THR > 0.01 and abs(sv - THR) / THR > 0.01
        agree = (pv < THR) == x["certifies"]
        pa.append(dict(family=x["family"], w=x["w"], k=x["k"], certificate=x["certificate"],
                       probe=f"{pv:.4e}", canonical=f"{sv:.4e}", probe_certifies=pv < THR,
                       canonical_certifies=x["certifies"], clear_margin=clear, agree=agree))
        if clear and not agree:
            problems.append(f"probe disagreement {key}")
    write(EXP / "results" / "probe_agreement.csv", pa)

    # ---- thresholds (verification step 6)
    tab_beyond = {}
    for x in rd(ROWS14):
        f, k = x["family"], int(x["k"])
        if k <= rc3.K_TAB_FIRST[f]:
            continue
        for c, src in list(REF14.items()) + [("exact", None)]:
            if src is not None and x["certificate"] != src:
                continue
            if src is None and x["certificate"] != "N1":
                continue
            lg = float(x["log10_S"] if src is None else x["log10_sigma"])
            if abs(lg - math.log10(THR)) < 1e-5:
                problems.append(f"tabulated row {f} k={k} {c} too close to the threshold")
            tab_beyond.setdefault((f, c), []).append((k, lg < math.log10(THR)))
    th = []
    for (f, w, c), v in sorted(verdicts.items(), key=lambda kv: (kv[0][0], Fraction(kv[0][1]), ALL.index(kv[0][2]))):
        reg, qt = meta[(f, w, min(v))]
        beyond = tab_beyond.get((f, c), []) if reg == "published" else []
        beyond_ok = all(ok for _, ok in beyond)
        complete = sorted(v) == sorted(rc3.K_PUBLISHED[f] if reg == "published" else rc3.K_EXT)
        kf, ks, fails = rc3.thresholds(v, beyond_ok)
        row = next(x for x in agg if (x["family"], x["w_exact"], x["certificate"]) == (f, w, c))
        th.append(dict(family=f, regime=reg, w=wdec(w), w_exact=w, q_target=row["q_target"],
                       q=row["q"], eta=row["eta"], mu=row["mu"], certificate=c,
                       k_scanned=f"{min(v)}-{max(v)}", scan_complete=complete,
                       k_first=kf if kf is not None else "",
                       k_star=ks if ks is not None else f">{max(v)}",
                       fails=";".join(map(str, fails)),
                       tab_rows_beyond=";".join(f"{k}:{'ok' if ok else 'FAIL'}" for k, ok in beyond)))
    onset = {}
    for f in rc3.EXT_FAMS:
        for q in rc3.Q_GRID:
            w = str(rc3.grid_w(f, q))
            ct, co = verdicts.get((f, w, "CT"), {}), verdicts.get((f, w, "contour"), {})
            hit = [k for k in ct if not ct[k] and co.get(k)]
            if hit and f not in onset:
                onset[f] = (q, w, hit)
    for x in th:
        if x["regime"] == "stronger" and x["family"] in onset:
            x["onset"] = x["w_exact"] == onset[x["family"]][1]
            x["strong"] = x["q_target"] == f"{float(rc3.Q_STRONG):.2f}"
    write(EXP / "results" / "thresholds.csv", th)

    # ---- dashboard contract
    figs = [dict(title="P4 modules, published parameters: log10 Sigma_1 against k "
                       "(certified iff below log10(sqrt2-1) = -0.383)",
                 plot_type="line", x_axis="k", y_axis="log10_sigma", z_axis=None,
                 series="certificate", filters={"family": "p4", "regime": "published"},
                 log_scale_x=False, log_scale_y=False)]
    for f, lab in (("cliques", "K4"), ("dimer", "dimer")):
        figs.append(dict(title=f"{lab} modules, stronger bridges q = 0.40: log10 Sigma_1 against k",
                         plot_type="line", x_axis="k", y_axis="log10_sigma", z_axis=None,
                         series="certificate", filters={"family": f, "q_target": "0.40"},
                         log_scale_x=False, log_scale_y=False))
    (EXP / "results" / "viz_schema.json").write_text(json.dumps(dict(figures=figs), indent=2) + "\n")

    # ---- report
    print(f"rows complete: {len(R)}/{len(meta)}")
    for x in rep:
        print(f"  repro {x['family']:7s} k={x['k']:2d} {x['certificate']:8s} diff {x['diff']:>13s}"
              f" {'ok' if x['passes'] else 'FAIL'}")
    for x in th:
        print(f"  {x['family']:7s} {x['regime']:9s} w={x['w']:>5s} q={x['q']} {x['certificate']:8s}"
              f" k_first={x['k_first']!s:>2} k*={x['k_star']!s:>4} fails=[{x['fails']}]"
              f" {x.get('tab_rows_beyond', '')}")
    print("onset:", {f: (str(q), w, hit) for f, (q, w, hit) in onset.items()})
    print(f"probe agreement: {sum(x['agree'] for x in pa)}/{len(pa)} "
          f"(clear-margin disagreements: {sum(x['clear_margin'] and not x['agree'] for x in pa)})")
    print("problems:", problems if problems else "none")
    if a.final and not problems:
        write_summary(agg, th, rep, pa, onset)
    return 0 if not problems else 1


def write_summary(agg, th, rep, pa, onset):
    p = EXP / "experiment_log.json"
    log = json.loads(p.read_text())
    runs = log["runs"]
    get = {(x["family"], x["regime"], x["w"], x["certificate"]): x for x in th}
    def kstar(f, reg, w):
        return {c: get[(f, reg, w, c)]["k_star"] for c in ALL}
    p4 = kstar("p4", "published", "0.3")
    groups = {}
    for r in runs:
        groups.setdefault(r["parameters"].get("tag", "?"), []).append(r["run_id"])
    ctr = rd(EXP / "results" / "ct_refine.csv") if (EXP / "results" / "ct_refine.csv").exists() else []
    log["summary"] = dict(
        n_runs_total=len(runs), n_runs_success=sum(r["status"] == "success" for r in runs),
        n_runs_failure=sum(r["status"] == "failure" for r in runs),
        run_groups={t: (ids[0] if len(ids) == 1 else f"{ids[0]}-{ids[-1]}") + f" ({len(ids)})"
                    for t, ids in groups.items()},
        ct_refine_check=dict(rows=len(ctr), all_still_fail=all(x["still_fails"] == "True" for x in ctr),
                             eta_capped_rows=sum(x["eta_capped_every_kappa"] == "True" for x in ctr)),
        aggregate_metrics=dict(
            rows=len({(x['family'], x['w'], x['k']) for x in agg}),
            violations=sum(int(x["violations"]) for x in agg),
            reproduction_pass=f"{sum(x['passes'] for x in rep)}/{len(rep)}",
            probe_clear_margin_disagreements=sum(x["clear_margin"] and not x["agree"] for x in pa),
            k_star_p4_published=p4,
            onset={f: dict(q=str(float(q)), w=w, ks_where_CT_fails_and_contour_certifies=hit)
                   for f, (q, w, hit) in onset.items()},
            k_star_q040={f: kstar(f, "stronger", f"{float(rc3.grid_w(f, rc3.Q_STRONG)):g}")
                         for f in rc3.EXT_FAMS}),
        matches_hypothesis="see README",
    )
    p.write_text(json.dumps(log, indent=2, ensure_ascii=False) + "\n")


if __name__ == "__main__":
    sys.exit(main())
