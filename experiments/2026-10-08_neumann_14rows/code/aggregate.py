"""Aggregate raw_<tag>.json into results/aggregate.csv and results/blocks.csv, with the
validity checks of plan.md (verification steps 4-6).

  python aggregate.py [--tag final] [--oracle]     # --oracle (re)builds probe_oracle.csv

Per-block bounds are exact rational upper endpoints, so every Sigma below is exact
rational arithmetic on certified upper bounds; mins over certificates are valid bounds.
Oracle (reproduction only, never a certificate):
  - the float probe of 2026-10-08 re-run at full digits (sondagem_n123.bounds, imported
    read-only) for N1, N2 and N3 (= min(N3, N2) in the probe);
  - the float N3x midpoint quadrature (1200 points per contour, as in sondagem_portas);
  - the probe's printed one-decimal tables, parsed, as a consistency check of the above.
"""
import argparse
import csv
import json
import math
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
EXP = HERE.parent
sys.path.insert(0, str(HERE))
import flint                                       # noqa: E402
from flint import arb, fmpq                        # noqa: E402
import neumann as nm                               # noqa: E402

ROOT = nm.ROOT
BENCH = ROOT / "experiments/2026-09-18_exact_benchmarks/results"
CERT = ROOT / "experiments/2026-09-18_interval_certification/results/aggregate.csv"
PROBE = ROOT / "derivacoes/2026-10-08_neumann_pontes"
CLOSED = ["N1", "N2", "N2p_stated", "N2p_path", "N2p_ladder", "N3", "N3vN2"]
ORDER = ["Sigma1"] + CLOSED + ["N3x_M1000", "N3x_M2000", "N3x", "best_closed", "best"]


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
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)


def lg(q):
    return nm.lg10(arb(q)) if q != 0 else -math.inf


# ------------------------------------------------------------------ oracle
def parse_printed():
    out = {}
    for line in (PROBE / "sondagem_n123_saida.txt").read_text().splitlines()[1:]:
        m = re.match(r"(\w+) (\d+)\s.*\|\s+(.*)\|\s+(\d+)", line)
        if m:
            v = [float(x) for x in m.group(3).split()]
            out[(m.group(1), int(m.group(2)))] = dict(S=v[0], floor=v[1], Sig1=v[2], N1=v[3],
                                                      N2=v[4], N3=v[5])
    for line in (PROBE / "sondagem_portas_saida.txt").read_text().splitlines():
        m = re.match(r"(\w+) (\d+).*Nc=\s*(\S+) N3x=\s*(\S+)", line)
        if m:
            out[(m.group(1), int(m.group(2)))].update(Nc=float(m.group(3)), N3x=float(m.group(4)))
    return out


def build_oracle(rows):
    import numpy as np
    prec = flint.ctx.prec
    sys.path.insert(0, str(PROBE))
    import sondagem_n123 as sp                     # sets ctx.prec = 256 on import
    from sondagem_portas import stadium
    printed = parse_printed()
    out = []
    for name, k in rows:
        flint.ctx.prec = 256
        P = sp.setup(name, k)
        o, *_ = sp.bounds(P)
        s, t = P["s"], P["t"]
        rec = dict(family=name, k=k)
        for key in ("N1", "N2", "N3"):
            sg = sum((o[key][(i, s)] * o[key][(i, t)] for i in range(k)), arb(0))
            rec[f"probe_{key}"] = f"{nm.lg10(sg):.6f}"
        if (name, k) in printed and "N3x" in printed[(name, k)]:
            g, delta = P["g"], P["delta"]
            eta = delta / (g - delta)
            sq = math.sqrt(1 - eta ** 2)
            absV = np.abs(P["Hout"])
            L = {}
            for i in range(k):
                zs, dl = stadium(*P["J"][i], g / 2, 1200)
                acc = {s: 0.0, t: 0.0}
                for z in zs:
                    Rabs = np.zeros((P["n"], P["n"]))
                    for V in P["parts"]:
                        Hj = P["H"][np.ix_(V, V)]
                        Rabs[np.ix_(V, V)] = np.abs(np.linalg.inv(Hj - z * np.eye(len(V))))
                    y = np.linalg.solve(np.eye(P["n"]) - absV @ Rabs, np.eye(P["n"])[:, [s, t]])
                    m = Rabs @ y
                    for c, x in enumerate((s, t)):
                        acc[x] += np.linalg.norm(m[P["parts"][i], c]) * dl
                for x in (s, t):
                    L[(i, x)] = 1.0 if P["comm"][x] == i else min(eta, acc[x] / (2 * math.pi) / sq)
            rec["probe_N3x_float1200"] = f"{math.log10(sum(L[(i, s)] * L[(i, t)] for i in range(k))):.6f}"
        for key, v in printed.get((name, k), {}).items():
            rec[f"printed_{key}"] = v
        out.append(rec)
        print(f"oracle {name} {k} done", flush=True)
    flint.ctx.prec = prec
    write(EXP / "results" / "probe_oracle.csv", out)


# ------------------------------------------------------------------ aggregate
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="final")
    ap.add_argument("--oracle", action="store_true")
    a = ap.parse_args()
    flint.ctx.prec = 128
    raw = json.loads((EXP / "results" / f"raw_{a.tag}.json").read_text())
    rows = [tuple(r) for r in raw["rows"]]
    if a.oracle or not (EXP / "results" / "probe_oracle.csv").exists():
        build_oracle(rows)
    oracle = {(r["family"], int(r["k"])): r for r in rd(EXP / "results" / "probe_oracle.csv")}
    # public export: final.csv is byte-identical to the published aggregate.csv
    final = {(r["family"], int(r["k"])): r for r in rd(BENCH / "aggregate.csv")}
    masses = {(r["family"], int(r["k"]), int(r["block"])): r for r in rd(BENCH / "masses.csv")}
    sig1 = {(r["family"], int(r["k"])): r for r in rd(CERT)}
    L = {}                                          # (name, k) -> cert -> (i, x) -> fmpq
    meta = {}
    for t in raw["tasks"]:
        if t["status"] != "success":
            continue
        s = t["spec"]
        key = (s[1], s[2])
        d = L.setdefault(key, {})
        if s[0] == "closed":
            for c, blk in t["res"].items():
                d[c] = {tuple(int(v) for v in kk.split(",")): fmpq(*[int(p) for p in val.split("/")])
                        if "/" in val else fmpq(int(val)) for kk, val in blk.items()}
            meta[key] = t["extra"]
        else:
            i, M = s[3], s[4]
            for x, val in t["res"].items():
                d.setdefault(f"N3x_M{M}", {})[(i, int(x))] = \
                    fmpq(*[int(p) for p in val.split("/")]) if "/" in val else fmpq(int(val))
    agg, blocks, problems = [], [], []
    for name, k in rows:
        if (name, k) not in meta:
            problems.append(f"{name} {k}: closed forms missing")
            continue
        d = L[(name, k)]
        P = nm.rig_setup(name, k)
        s, t = P["s"], P["t"]
        keys = [(i, x) for i in range(k) for x in (s, t)]
        n3 = [c for c in ("N3x_M1000", "N3x_M2000") if c in d and len(d[c]) == len(keys)]
        if n3:
            d["N3x"] = {kk: min(d[c][kk] for c in n3) for kk in keys}
        d["best_closed"] = {kk: min(d[c][kk] for c in CLOSED if c in d) for kk in keys}
        d["best"] = {kk: min(d[c][kk] for c in CLOSED + n3 if c in d) for kk in keys}
        S = arb(final[(name, k)]["S"])
        mass = {(i, x): arb(masses[(name, k, i)]["mass_" + ("s" if x == s else "t")])
                for (i, x) in keys}
        floor = sum((mass[(i, s)] * mass[(i, t)] for i in range(k)), arb(0))
        tol = fmpq(999999, 1000000)
        orc = oracle.get((name, k), {})
        for c in ORDER:
            if c == "Sigma1":
                sg = fmpq(*[int(p) for p in str(nm.upper_q(arb(sig1[(name, k)]["Sigma1_upper_certified"]))).split("/")])
                viol, per_block = (1 if sg < nm.lower_q(S) else 0), None
            elif c in d:
                per_block = d[c]
                sg = sum((per_block[(i, s)] * per_block[(i, t)] for i in range(k)), fmpq(0))
                viol = sum(per_block[kk] < nm.lower_q(mass[kk]) * tol for kk in keys)
                viol += int(sg < nm.lower_q(S) * tol)
            else:
                continue
            rec = dict(family=name, k=k, D=meta[(name, k)]["D"], certificate=c,
                       log10_sigma=f"{lg(sg):.6f}", log10_S=f"{nm.lg10(S):.6f}",
                       log10_floor=f"{nm.lg10(floor):.6f}",
                       dist_S=f"{lg(sg) - nm.lg10(S):.6f}",
                       dist_floor=f"{lg(sg) - nm.lg10(floor):.6f}", violations=viol)
            pk = {"N1": "probe_N1", "N2": "probe_N2", "N3vN2": "probe_N3",
                  "N3x_M1000": "probe_N3x_float1200", "N3x_M2000": "probe_N3x_float1200"}.get(c)
            if pk and orc.get(pk):
                rec["probe_log10"] = orc[pk]
                rec["probe_diff"] = f"{lg(sg) - float(orc[pk]):.6f}"
            pp = {"N1": "printed_N1", "N2": "printed_N2", "N3vN2": "printed_N3",
                  "N2p_ladder": "printed_Nc", "N3x": "printed_N3x", "Sigma1": "printed_Sig1"}.get(c)
            if pp and orc.get(pp):
                rec["printed_log10"] = orc[pp]
            if c in ("N2", "N3", "N3vN2"):
                rec["rho_certified"] = meta[(name, k)]["cert"]["N2" if c == "N2" else "N3"]
            if c.startswith("N3x_M"):
                rec["rho_certified"] = True             # NotMMatrix would have failed the task
            if c == "N2p_ladder":
                rec["ladder_certified"] = meta[(name, k)]["info"]["ladder_certified"]
            agg.append(rec)
            if viol:
                problems.append(f"{name} {k} {c}: {viol} violations")
            if per_block is not None:
                for i in range(k):
                    blocks.append(dict(family=name, k=k, block=i, certificate=c,
                                       log10_L_s=f"{lg(per_block[(i, s)]):.6f}",
                                       log10_L_t=f"{lg(per_block[(i, t)]):.6f}",
                                       log10_mass_s=f"{nm.lg10(mass[(i, s)]):.6f}",
                                       log10_mass_t=f"{nm.lg10(mass[(i, t)]):.6f}"))
    tag = "" if a.tag == "final" else f"_{a.tag}"
    write(EXP / "results" / f"aggregate{tag}.csv", agg)
    write(EXP / "results" / f"blocks{tag}.csv", blocks)
    # ---- report
    print(f"{'row':11s} " + " ".join(f"{c:>9s}" for c in ORDER) + "      S")
    for name, k in rows:
        vals = {r["certificate"]: float(r["log10_sigma"]) for r in agg
                if (r["family"], r["k"]) == (name, k)}
        S = next(float(r["log10_S"]) for r in agg if (r["family"], r["k"]) == (name, k))
        print(f"{name+' '+str(k):11s} " + " ".join(f"{vals.get(c, float('nan')):9.2f}"
                                                  for c in ORDER) + f" {S:8.2f}")
    worst = max((abs(float(r["probe_diff"])) for r in agg
                 if r.get("probe_diff") and not r["certificate"].startswith("N3x")), default=0)
    print(f"max |closed-form - probe| = {worst:.2e} (tolerance 0.05)")
    n3d = [float(r["probe_diff"]) for r in agg if r.get("probe_diff") and r["certificate"] == "N3x_M2000"]
    if n3d:
        print(f"N3x(M=2000) - float N3x: min {min(n3d):+.4f}, max {max(n3d):+.4f} "
              f"(must be >= 0 up to quadrature error; tolerance +0.02)")
    print("problems:", problems if problems else "none")
    return 0 if not problems else 1


if __name__ == "__main__":
    sys.exit(main())
