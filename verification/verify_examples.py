"""Re-derives every number in the three tables of Section 9 of main.tex.

Claims tested: delta, g, W_i, eta, ghat, rho_i, d, mu, D_0, dist(s,V_i),
D=dist(s,t), Sigma_1, F_esp = min{1,Sigma_1^2}, B^2, and b_eta^2 --
plus a refutation search: the true S, Fbar and a densely sampled Fmax must
never exceed the certificate.
"""
import numpy as np
from build_from_scratch import (assemble, family_cliques, family_dimer,
                                family_p4, separation_parameters, certificate,
                                S_overlap, Fbar_exact, amplitude)

np.set_printoptions(precision=12, suppress=False)


def sup4(x):
    """4 significant figures, rounded UP (the paper presents upper bounds)."""
    if x == 0:
        return 0.0
    import math
    e = math.floor(math.log10(abs(x)))
    f = 10 ** (e - 3)
    return math.ceil(x / f) * f


def report(name, fam_fn, rows, s_fn, t_fn, paper_par, paper_rows):
    print("=" * 78)
    print(name)
    print("=" * 78)
    ok = True
    # ---- structural parameters, from the k of the first row
    k0 = rows[0][0]
    n, parts, comm_of, H, Dmat, d, mu = assemble(fam_fn(k0))
    par = separation_parameters(H, parts, comm_of, d, mu)
    D0 = max(max(Dmat[u, v] for u in Vi for v in Vi) for Vi in parts)
    print(f"  measured: delta={par['delta']:.10f}  g={par['g']:.10f}  "
          f"W_i={par['W'][0]:.10f}")
    print(f"            eta={par['eta']:.10f}  ghat={par['ghat']:.10f}  "
          f"rho_i={par['rho'][0]:.10f}  D_0={D0:.0f}  d={d} mu={mu}")
    print(f"  paper   : {paper_par}")
    for key, claimed in paper_par.items():
        got = {'delta': par['delta'], 'g': par['g'], 'W': par['W'][0],
               'eta': par['eta'], 'ghat': par['ghat'], 'rho': par['rho'][0],
               'D0': D0}[key]
        agree = abs(got - claimed) <= 5e-7 * max(1.0, abs(claimed))
        ok &= agree
        print(f"    {key:6s} paper={claimed!r:>14}  mine={got:.10f}  "
              f"{'OK' if agree else '*** MISMATCH ***'}")
    # b_eta^2, k-independent
    cert0 = certificate(H, parts, comm_of, d, mu, s_fn(k0), t_fn(k0), 0.5,
                        Dmat, par)
    print(f"  b_eta^2 = {cert0['b_eta']**2:.10e}  (4sf up: {sup4(cert0['b_eta']**2):.4e})")

    # ---- per-row table reproduction
    print(f"  {'k':>3} {'n':>4} {'D':>4} {'kappa':>7} | {'Sigma1 mine':>13} "
          f"{'paper':>11} | {'F_esp mine':>12} {'paper':>11} | "
          f"{'B^2 mine':>12} {'paper':>11} | true S    bound-ok")
    for (k, kappa, pS1, pFe, pB2) in paper_rows:
        n, parts, comm_of, H, Dmat, d, mu = assemble(fam_fn(k))
        s, t = s_fn(k), t_fn(k)
        c = certificate(H, parts, comm_of, d, mu, s, t, kappa, Dmat)
        # distances claimed in the text
        dist_sVi = [min(Dmat[s, u] for u in Vi) for Vi in parts]
        dist_tVi = [min(Dmat[t, u] for u in Vi) for Vi in parts]
        step = {'cliques': 2, 'dimer': 2, 'p4': 4}[name.split()[0]]
        dist_ok = (all(dist_sVi[i] == step * i for i in range(k)) and
                   all(dist_tVi[i] == step * (k - 1 - i) for i in range(k)))
        Dclaim = step * k - 1
        dist_ok &= (c['D'] == Dclaim)
        trueS = S_overlap(H, s, t)
        trueFbar = Fbar_exact(H, s, t)
        # dense refutation sweep for Fmax
        taus = np.linspace(0, 400.0, 400001)
        Fs = np.abs(amplitude(H, s, t, taus)) ** 2
        trueFmax_lb = Fs.max()
        m1 = abs(sup4(c['Sigma1']) - pS1) <= 1.0001e-4 * abs(pS1)
        m2 = abs(sup4(c['F_esp']) - pFe) <= 1.0001e-4 * abs(pFe)
        m3 = abs(sup4(c['B2']) - pB2) <= 1.0001e-4 * abs(pB2)
        bound_ok = (trueS <= c['B'] + 1e-12 and trueFmax_lb <= c['B2'] + 1e-12
                    and trueFbar <= min(c['Sigma2'], 0.5 * c['B'] ** 2) + 1e-12)
        ok &= m1 and m2 and m3 and dist_ok and bound_ok
        print(f"  {k:>3} {n:>4} {int(c['D']):>4} {kappa:>7.4f} | "
              f"{sup4(c['Sigma1']):>13.4e} {pS1:>11.4e} {'ok' if m1 else 'XX'} | "
              f"{sup4(c['F_esp']):>12.4e} {pFe:>11.4e} {'ok' if m2 else 'XX'} | "
              f"{sup4(c['B2']):>12.4e} {pB2:>11.4e} {'ok' if m3 else 'XX'} | "
              f"{trueS:.3e} {'dist-ok' if dist_ok else 'DIST-XX'} "
              f"{'bnd-ok' if bound_ok else 'BOUND-VIOLATED'}")
    print(f"  --> {name}: {'ALL REPRODUCED' if ok else 'DISCREPANCY'}")
    return ok


all_ok = True
all_ok &= report(
    "cliques  (Example 9.1, K_4 modules)", family_cliques,
    [(3,), (4,), (5,), (6,), (8,), (10,)],
    lambda k: 0, lambda k: 4 * k - 1,
    dict(delta=0.10, g=20.0, W=4.0, eta=0.0050251, ghat=2.4750,
         rho=1.135041, D0=1),
    [(3, 0.5000, 0.01008, 1.016e-4, 1.011e-4),
     (4, 0.5000, 0.01011, 1.021e-4, 1.011e-4),
     (5, 0.8440, 3.588e-3, 1.288e-5, 1.288e-5),
     (6, 0.8712, 4.279e-4, 1.831e-7, 1.831e-7),
     (8, 0.8615, 2.353e-5, 5.535e-10, 5.535e-10),
     (10, 0.8681, 8.805e-7, 7.752e-13, 7.752e-13)])

all_ok &= report(
    "dimer  (Example 9.2, detuned dimerized chain)", family_dimer,
    [(3,), (5,), (7,), (9,), (12,)],
    lambda k: 0, lambda k: 2 * k - 1,
    dict(delta=0.05, g=24 - 2 * np.sqrt(26), W=2 * np.sqrt(26),
         eta=0.0036358, ghat=3.425490, rho=1.478467, D0=1),
    [(3, 0.5000, 7.285e-3, 5.307e-5, 5.288e-5),
     (5, 0.8565, 5.477e-4, 3.000e-7, 3.000e-7),
     (7, 0.8830, 1.721e-5, 2.960e-10, 2.960e-10),
     (9, 0.8739, 1.983e-7, 3.931e-14, 3.931e-14),
     (12, 0.9006, 9.191e-11, 8.448e-21, 8.448e-21)])

all_ok &= report(
    "p4  (Example 9.3, P_4 modules)", family_p4,
    [(10,), (20,), (40,)],
    lambda k: 0, lambda k: 4 * k - 1,
    dict(delta=0.30, g=4 - np.sqrt(5), W=1 + np.sqrt(5),
         eta=0.2049275, ghat=0.290983, rho=3.098161, D0=3),
    [(10, 0.8642, 0.2474, 0.06118, 0.06118),
     (20, 0.8951, 7.449e-4, 5.548e-7, 5.548e-7),
     (40, 0.9453, 1.290e-11, 1.662e-22, 1.662e-22)])

print("\nOVERALL:", "ALL EXAMPLE CLAIMS REPRODUCED" if all_ok else "SEE DISCREPANCIES ABOVE")
