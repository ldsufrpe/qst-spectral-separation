"""Fits, predictions, rate shares and end effect for the Stark test (plan.md, hypothesis).

Reads results/z_series.csv and results/chains.csv; writes results/fits.csv,
results/rate_shares.csv, results/chain_summary.csv and the long-format
results/aggregate.csv (the dashboard contract).

Predictions (fixed before the canonical run; the 'stated' one after the toy):
  per-bridge N2'   c = log10 B(q'^2),     q' = 2J/F,  b = 1/2
  ladder N2'       c = log10(2 B(q^2)),   q  = 4J/F,  b = 1/2
  stated N2'       c = log10(4 B(q^2))    (derived from the toy's F/J = 10 value, not before it)
Fit: log10(bound/mass) = a + b log10 N + c N over N = 10..100 (eta floor excluded).
"""
import csv
import math
from pathlib import Path

import numpy as np

EXP = Path(__file__).resolve().parent.parent
FORMS = ("Sigma1", "N1", "N2p_stated", "N2p_ladder", "N2p_bridge")


def B(z):
    return (1 - math.sqrt(1 - 4 * z)) / (2 * z)


def frac(s):
    if "/" in s:
        p, q = s.split("/")
        return int(p) / int(q)
    return float(s)


def write(path, rows):
    keys = []
    for r in rows:
        keys += [k for k in r if k not in keys]
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)


def predictions(F):
    q, qb = 4 / F, 2 / F
    out = {"N2p_bridge": (math.log10(B(qb ** 2)), 0.5) if 4 * qb ** 2 < 1 else None}
    if 4 * q ** 2 < 1:
        out["N2p_ladder"] = (math.log10(2 * B(q ** 2)), 0.5)
        out["N2p_stated"] = (math.log10(4 * B(q ** 2)), None)
    return out


def main():
    z = list(csv.DictReader(open(EXP / "results" / "z_series.csv")))
    Fs = list(dict.fromkeys(r["F_over_J"] for r in z))
    fits, shares, agg = [], [], []
    for Fs_ in Fs:
        F = frac(Fs_)
        eta = 2 / (F - 2)
        rows = [r for r in z if r["F_over_J"] == Fs_]
        pred = predictions(F)
        for r in rows:
            for form in FORMS:
                if r.get(f"{form}_log10_Sigma"):
                    agg.append(dict(part="Z", F_over_J=Fs_, N=int(r["N"]), certificate=form,
                                    log10_sigma=r[f"{form}_log10_Sigma"], log10_S=r["log10_S"],
                                    dist_S=f"{float(r[f'{form}_log10_Sigma']) - float(r['log10_S']):.6f}",
                                    log10_mass=r["log10_mass"], log10_L=r[f"{form}_log10_L"]))
        for form in FORMS:
            pts = [(int(r["N"]), float(r[f"{form}_log10_L"]) - float(r["log10_mass"]))
                   for r in rows if r.get(f"{form}_log10_L") and 10 <= int(r["N"]) <= 100
                   and float(r[f"{form}_log10_L"]) < math.log10(eta) - 1e-9]
            if len(pts) < 5:
                continue
            N = np.array([p[0] for p in pts], float)
            y = np.array([p[1] for p in pts])
            A = np.vstack([np.ones_like(N), np.log10(N), N]).T
            coef, *_ = np.linalg.lstsq(A, y, rcond=None)
            rms = float(np.sqrt(np.mean((A @ coef - y) ** 2)))
            A2 = np.vstack([np.ones_like(N), np.log10(N), N, 1 / N]).T
            coef2, *_ = np.linalg.lstsq(A2, y, rcond=None)
            rec = dict(F_over_J=Fs_, form=form, n_points=len(pts), a=f"{coef[0]:.6f}",
                       b=f"{coef[1]:.6f}", c=f"{coef[2]:.8f}", rms=f"{rms:.3e}",
                       c_with_1overN=f"{coef2[2]:.8f}", b_with_1overN=f"{coef2[1]:.6f}")
            if pred.get(form):
                cp, bp = pred[form]
                rec.update(c_pred=f"{cp:.8f}", c_absdev=f"{coef[2] - cp:+.3e}",
                           c_reldev=f"{(coef[2] - cp) / cp:+.4f}",
                           c_pass_2pct=abs(coef[2] - cp) <= 0.02 * abs(cp))
                if bp is not None:
                    rec.update(b_pred=bp, b_pass_005=abs(coef[1] - bp) <= 0.05)
            fits.append(rec)
        # S-level law for the per-bridge form
        if pred.get("N2p_bridge"):
            cp = pred["N2p_bridge"][0]
            dev = [float(r["N2p_bridge_log10_Sigma"]) - float(r["log10_S"])
                   - (math.log10(math.pi * int(r["N"]) / 2) + int(r["N"]) * cp)
                   for r in rows if 10 <= int(r["N"]) <= 100]
            fits.append(dict(F_over_J=Fs_, form="N2p_bridge_S_level", n_points=len(dev),
                             a=f"{np.mean(dev):.6f}", rms=f"{np.std(dev):.3e}",
                             spread=f"{max(dev) - min(dev):.4f}",
                             spread_pass_05=(max(dev) - min(dev)) <= 0.5))
        # rate shares over N-windows
        byN = {int(r["N"]): r for r in rows}
        for (n1, n2) in ((10, 20), (40, 50), (90, 100)):
            dS = float(byN[n2]["log10_S"]) - float(byN[n1]["log10_S"])
            for form in FORMS:
                k = f"{form}_log10_Sigma"
                if byN[n1].get(k) and byN[n2].get(k):
                    shares.append(dict(part="Z", F_over_J=Fs_, window=f"{n1}-{n2}", certificate=form,
                                       share=f"{(float(byN[n2][k]) - float(byN[n1][k])) / dS:.4f}"))
    # ---- chains
    summ = []
    cpath = EXP / "results" / "chains.csv"
    if cpath.exists():
        ch = list(csv.DictReader(open(cpath)))
        certs = ("Sigma1", "N1", "N2", "N2p_stated", "N2p_ladder", "N2p_path", "N3x")
        for r in ch:
            F, N = frac(r["F_over_J"]), int(r["N"])
            cp = math.log10(B((2 / F) ** 2))
            pred = math.log10(math.pi * N / 2) + N * cp
            rec = dict(L=r["L"], N=N, F_over_J=r["F_over_J"], log10_S=r["log10_S"],
                       end_effect=f"{float(r['log10_S']) - float(r['log10_S_Z']):.6f}",
                       pred_bridge_S_level=f"{pred:.6f}", dps=r["dps"],
                       violations_total=r["violations_total"], N2_le_N2p_path=r["N2_le_N2p_path"])
            for c in certs:
                if r.get(f"{c}_log10_Sigma"):
                    d = float(r[f"{c}_log10_Sigma"]) - float(r["log10_S"])
                    rec[f"dist_{c}"] = f"{d:.6f}"
                    agg.append(dict(part="chain", F_over_J=r["F_over_J"], N=N, certificate=c,
                                    log10_sigma=r[f"{c}_log10_Sigma"], log10_S=r["log10_S"],
                                    dist_S=f"{d:.6f}"))
            rec["N2_within_pred_plus1"] = float(rec["dist_N2"]) <= pred + 1
            summ.append(rec)
        for Fs_ in dict.fromkeys(r["F_over_J"] for r in ch):
            byL = {int(r["L"]): r for r in ch if r["F_over_J"] == Fs_}
            for (l1, l2) in ((21, 61), (41, 61)):
                if l1 in byL and l2 in byL:
                    dS = float(byL[l2]["log10_S"]) - float(byL[l1]["log10_S"])
                    for c in certs:
                        k = f"{c}_log10_Sigma"
                        if byL[l1].get(k) and byL[l2].get(k):
                            shares.append(dict(part="chain", F_over_J=Fs_, window=f"L{l1}-L{l2}",
                                               certificate=c,
                                               share=f"{(float(byL[l2][k]) - float(byL[l1][k])) / dS:.4f}"))
    write(EXP / "results" / "fits.csv", fits)
    write(EXP / "results" / "rate_shares.csv", shares)
    if summ:
        write(EXP / "results" / "chain_summary.csv", summ)
    write(EXP / "results" / "aggregate.csv", agg)
    # ---- report
    print("Z fits (N=10..100): form, F/J: b, c, c_pred, reldev, pass")
    for f in fits:
        if f["form"] in ("N2p_bridge", "N2p_ladder", "N2p_stated"):
            print(f"  {f['form']:11s} {f['F_over_J']:>5s}: b={float(f['b']):+.3f} c={float(f['c']):.6f} "
                  f"pred={f.get('c_pred', '-'):>10s} rel={f.get('c_reldev', '-'):>8s} "
                  f"abs={f.get('c_absdev', '-'):>10s} c2pct={f.get('c_pass_2pct', '-')} "
                  f"b={f.get('b_pass_005', '-')} | with 1/N: c={float(f['c_with_1overN']):.6f}")
        elif f["form"] == "N2p_bridge_S_level":
            print(f"  S-level bridge {f['F_over_J']:>5s}: offset {float(f['a']):+.3f}, spread "
                  f"{f['spread']} pass={f['spread_pass_05']}")
    print("rate shares:")
    for s in shares:
        if s["certificate"] in ("N2p_bridge", "N2", "N3x", "N1", "Sigma1", "N2p_path"):
            print(f"  {s['part']:5s} F/J={s['F_over_J']:>5s} {s['window']:>8s} {s['certificate']:11s} {s['share']}")
    for r in summ:
        print(f"  chain L={r['L']:>2s} F/J={r['F_over_J']:>5s}: end {float(r['end_effect']):+.3f} "
              f"dist N2 {float(r['dist_N2']):6.2f} (pred {float(r['pred_bridge_S_level']):6.2f}) "
              f"N3x {r.get('dist_N3x', '—'):>9s} N1 {float(r['dist_N1']):6.1f} "
              f"Sig1 {float(r['dist_Sigma1']):6.1f} viol {r['violations_total']}")


if __name__ == "__main__":
    main()
