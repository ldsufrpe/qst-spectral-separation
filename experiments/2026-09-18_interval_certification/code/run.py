"""Certify every entry of Tables Tab:cliques, Tab:dimero, Tab:p4 and the displayed
parameters of Examples 7.1-7.3 in Arb ball arithmetic.

Usage:  python run.py [--smoke] [--prec BITS]
The tables and displayed parameters are parsed at run time from data/examples_sec7.tex,
a verbatim extract of Examples 7.1-7.3 (tables included) of the manuscript.
"""
import argparse
import csv
import re
import sys
import time
from pathlib import Path

import flint
from flint import arb, fmpq

import certify as C
import logutil as LU
from families import build

MAIN = LU.EXP / "data" / "examples_sec7.tex"
RES = LU.EXP / "results"
PROV = ("Addresses delegated review finding R-C4 (revision_response_trace.json): "
        "reproducible rigorous certification of every table entry of Section 7, "
        "including the decimal kappa selection; library/version/precision named here.")

TABLES = {"cliques": "Tab:cliques", "dimer": "Tab:dimero", "p4": "Tab:p4"}

# Parameters displayed in the running text of Examples 7.1-7.3 ("≈" = nearest rounding)
DISPLAYED = [
    ("cliques", "eta", "0.0050251"), ("cliques", "ghat", "2.4750"),
    ("cliques", "rho", "1.135041"), ("cliques", "rho", "1.135"),
    ("dimer", "W", "10.198039"), ("dimer", "g", "13.801961"), ("dimer", "g", "13.802"),
    ("dimer", "eta", "0.0036358"), ("dimer", "ghat", "3.425490"),
    ("dimer", "rho", "1.478467"),
    ("p4", "g", "1.763932"), ("p4", "eta", "0.2049275"), ("p4", "ghat", "0.290983"),
    ("p4", "rho", "3.098161"), ("p4", "rho", "3.10"), ("p4", "W", "3.24"),
    ("p4", "g-2delta", "1.16"),
]
# Displayed upper bounds b_eta^2 <= value (4 s.f.)
DISPLAYED_BOUNDS = [("cliques", "1.011\\times10^{-4}"), ("dimer", "5.288\\times10^{-5}"),
                    ("p4", "0.1610")]


def parse_num(tex: str) -> fmpq:
    tex = tex.strip().strip("$").replace(" ", "").replace("\\\\", "")
    m = re.fullmatch(r"([0-9.]+)(?:\\times10\^\{(-?\d+)\})?", tex)
    if not m:
        raise ValueError(tex)
    mant, e = m.group(1), int(m.group(2) or 0)
    dec = len(mant.split(".")[1]) if "." in mant else 0
    val = fmpq(int(mant.replace(".", "")), 10 ** dec)
    return val * C.pow10(e)


def parse_table(tex, label):
    body = tex.split("\\label{" + label + "}")[1].split("\\bottomrule")[0]
    rows = []
    for line in body.splitlines():
        if not re.match(r"^\$\d+\$ &", line.strip()):
            continue
        cells = [c.strip() for c in line.strip().rstrip("\\").split("&")]
        rows.append(dict(k=int(cells[0].strip("$")), n=int(cells[1].strip("$")),
                         D=int(cells[2].strip("$")), kappa_str=cells[3].strip("$"),
                         kappa=parse_num(cells[3]), Sigma1=cells[4], F_esp=cells[5],
                         B2=cells[6].rstrip("\\").strip()))
    return rows


def fmt(x: fmpq, digits=6):
    return f"{float(x.p) / float(x.q):.{digits}e}" if x != 0 else "0"


def kappa_audit(G, par, geo, kappa_tab, step_den):
    """Sigma1 on the decimal grid kappa = j/step_den, j=1..step_den-1.
    Returns the ball-aware argmin set and whether kappa_tab belongs to it."""
    his, los = [], []
    for j in range(1, step_den):
        b = C.sigma1(G, par, geo, fmpq(j, step_den))
        his.append(C.upper_q(b))
        los.append(C.lower_q(b))
    mhi = min(his)
    argset = [j for j in range(1, step_den) if los[j - 1] <= mhi]
    jtab = kappa_tab * step_den
    in_set = jtab.q == 1 and int(jtab.p) in argset
    return dict(argmin_first=fmpq(argset[0], step_den), argmin_last=fmpq(argset[-1], step_den),
                argmin_size=len(argset), grid_min_upper=mhi, kappa_tab_in_argmin=in_set,
                sigma1_at_tab_upper=his[int(jtab.p) - 1] if jtab.q == 1 else None)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--prec", type=int, default=256)
    ap.add_argument("--no-kappa", action="store_true")
    args = ap.parse_args()
    C.set_prec(args.prec)
    t0, ts = time.time(), LU.now()
    info = LU.capture_env()
    hdr = LU.header("2026-09-18_interval_certification", "symbolic", ["symbolic"], PROV,
                    "2026-09-18T00:00:00-03:00")
    tex = MAIN.read_text()
    RES.mkdir(exist_ok=True)
    tag = "smoke" if args.smoke else f"prec{args.prec}"
    step_den = 100 if args.smoke else 10000

    rows_out, audit_out, par_out, failures = [], [], [], []
    params_cache = {}
    for fam, label in TABLES.items():
        rows = parse_table(tex, label)
        if args.smoke:
            rows = rows[:1]
        for r in rows:
            G = build(fam, r["k"])
            assert G["n"] == r["n"], (fam, r)
            par = C.parameters(G)
            params_cache.setdefault(fam, (G, par))
            geo = C.geometry(G)
            assert geo["D"] == r["D"], (fam, r["k"], geo["D"], r["D"])
            cert = C.row_certificates(G, par, geo, r["kappa"])
            out = dict(family=fam, k=r["k"], n=r["n"], D=r["D"], kappa_tab=r["kappa_str"])
            for col in ("Sigma1", "F_esp", "B2"):
                printed = parse_num(r[col])
                ok_u, ok_t, ceil4 = C.check_printed(printed, cert[col])
                out[f"{col}_printed"] = fmt(printed, 3)
                out[f"{col}_upper_certified"] = fmt(C.upper_q(cert[col]), 12)
                out[f"{col}_ball_radius"] = fmt(C.upper_q(cert[col]) - C.lower_q(cert[col]), 2)
                out[f"{col}_ceil4"] = fmt(ceil4, 3)
                out[f"{col}_ok_upper"] = ok_u
                out[f"{col}_ok_tight"] = ok_t
                if not ok_u:
                    failures.append(f"{fam} k={r['k']} {col}: printed {fmt(printed,3)} < certified upper {fmt(C.upper_q(cert[col]),6)}")
            out["active_bound"] = ("Sigma1" if C.upper_q(cert["B2"]) == C.upper_q(cert["F_esp"])
                                   else "b_eta" if C.upper_q(cert["B2"]) == C.upper_q(cert["b_eta2"])
                                   else "B_E")
            if not args.no_kappa:
                ka = kappa_audit(G, par, geo, r["kappa"], step_den)
                out["kappa_tab_in_grid_argmin"] = ka["kappa_tab_in_argmin"]
                out["kappa_argmin_range"] = f"[{fmt(ka['argmin_first'],4)}, {fmt(ka['argmin_last'],4)}]"
                out["kappa_argmin_size"] = ka["argmin_size"]
                out["Sigma1_grid_min_upper"] = fmt(ka["grid_min_upper"], 6)
                audit_out.append(dict(family=fam, k=r["k"], kappa_tab=r["kappa_str"],
                                      in_argmin=ka["kappa_tab_in_argmin"],
                                      argmin_first=fmt(ka["argmin_first"], 4),
                                      argmin_last=fmt(ka["argmin_last"], 4),
                                      argmin_size=ka["argmin_size"],
                                      sigma1_grid_min_upper=fmt(ka["grid_min_upper"], 6),
                                      sigma1_at_kappa_tab_upper=fmt(ka["sigma1_at_tab_upper"], 6)
                                      if ka["sigma1_at_tab_upper"] is not None else ""))
            rows_out.append(out)
            print(f"{fam:8s} k={r['k']:3d} " + " ".join(
                f"{c}:{'OK' if out[c+'_ok_upper'] else 'FAIL'}/{'tight' if out[c+'_ok_tight'] else 'loose'}"
                for c in ("Sigma1", "F_esp", "B2")) +
                (f" kappa_in_argmin={out.get('kappa_tab_in_grid_argmin')}" if not args.no_kappa else ""),
                flush=True)

    # displayed parameters
    for fam, name, s in ([] if args.smoke else DISPLAYED):
        G, par = params_cache[fam]
        val = {"eta": par["eta"], "ghat": par["ghat"], "rho": par["rho"][0],
               "W": par["W"][0], "g": par["g"], "g-2delta": par["g"] - 2 * par["delta"]}[name]
        ok = C.rounds_to(s, val)
        present = s in tex
        par_out.append(dict(family=fam, parameter=name, printed=s, certified_mid=fmt(C.upper_q(val), 12),
                            rounds_correctly=ok, string_found_in_main_tex=present))
        if not ok:
            failures.append(f"displayed {fam} {name}={s} does not round from {fmt(C.upper_q(val),10)}")
    for fam, s in ([] if args.smoke else DISPLAYED_BOUNDS):
        G, par = params_cache[fam]
        printed = parse_num(s)
        ok_u, ok_t, c4 = C.check_printed(printed, C.b_eta_sq(par))
        par_out.append(dict(family=fam, parameter="b_eta^2 (upper bound)", printed=s,
                            certified_mid=fmt(C.upper_q(C.b_eta_sq(par)), 12),
                            rounds_correctly=ok_u and ok_t, string_found_in_main_tex=s in tex))
        if not ok_u:
            failures.append(f"b_eta^2 {fam}: printed {s} below certified upper")
    charpolys = {fam: par["charpolys"][:2] for fam, (G, par) in params_cache.items()}

    # write results
    agg = RES / ("aggregate.csv" if not args.smoke else "smoke_aggregate.csv")
    if tag == "prec512":
        agg = RES / "aggregate_prec512.csv"
    with open(agg, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows_out[0].keys()))
        w.writeheader()
        w.writerows(rows_out)
    outs = [str(agg.relative_to(LU.EXP))]
    if par_out:
        p = RES / (f"params_{tag}.csv" if tag != "prec256" else "params.csv")
        with open(p, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(par_out[0].keys()))
            w.writeheader()
            w.writerows(par_out)
        outs.append(str(p.relative_to(LU.EXP)))
    if audit_out and not args.smoke:
        p = RES / (f"kappa_audit_{tag}.csv" if tag != "prec256" else "kappa_audit.csv")
        with open(p, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(audit_out[0].keys()))
            w.writeheader()
            w.writerows(audit_out)
        outs.append(str(p.relative_to(LU.EXP)))

    n_ok = sum(all(r[f"{c}_ok_upper"] for c in ("Sigma1", "F_esp", "B2")) for r in rows_out)
    n_tight = sum(all(r[f"{c}_ok_tight"] for c in ("Sigma1", "F_esp", "B2")) for r in rows_out)
    status = "success" if not failures else "failure"
    run = dict(timestamp_start=ts, timestamp_end=LU.now(),
               wallclock_seconds=round(time.time() - t0, 1), peak_memory_mb=round(LU.peak_mb(), 1),
               status=status, seed=None,
               parameters=dict(scale="smoke" if args.smoke else "full", precision_bits=args.prec,
                               kappa_grid_step=f"1/{step_den}", library=f"python-flint {flint.__version__}",
                               manuscript=str(MAIN.relative_to(LU.ROOT)),
                               manuscript_sha256=LU.sha256_file(MAIN)),
               inputs=dict(data_files=[str(MAIN.relative_to(LU.ROOT))],
                           data_hashes={MAIN.name: LU.sha256_file(MAIN)}),
               outputs=dict(result_files=outs,
                            result_hashes={o: LU.sha256_file(LU.EXP / o) for o in outs},
                            scalar_results=dict(rows=len(rows_out), rows_all_upper_ok=n_ok,
                                                rows_all_tight=n_tight,
                                                kappa_in_argmin=sum(bool(r.get("kappa_tab_in_grid_argmin")) for r in rows_out),
                                                ball_overlaps_in_min=len(C.OVERLAPS),
                                                block_charpolys_sample=charpolys),
                            failures=failures),
               exit_code=0 if not failures else 1,
               provenance_note=PROV,
               notes="seed=null: deterministic exact/ball arithmetic, no randomness.")
    rid = LU.append_run(hdr, run)
    print(f"{rid}: {status}  rows upper-ok {n_ok}/{len(rows_out)}  tight {n_tight}/{len(rows_out)}  "
          f"overlaps {len(C.OVERLAPS)}  {run['wallclock_seconds']} s")
    for f_ in failures:
        print("  FAIL:", f_)
    return 0 if not failures else 1


if __name__ == "__main__":
    sys.exit(main())
