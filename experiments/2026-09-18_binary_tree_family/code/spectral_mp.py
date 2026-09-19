"""High-precision exact spectral quantities, used only to CALIBRATE the certificates
(never to produce them).

All eigen-decompositions are mpmath `eigsy` at a chosen number of decimal digits.
The certificates go down to 1e-22; float64 eigenvectors cannot resolve projector
entries below ~1e-16 ||H||, hence the multiprecision.
"""
import math

import mpmath as mp
import numpy as np
from scipy.optimize import minimize_scalar


def to_mp(G_dense):
    n = len(G_dense)
    H = mp.matrix(n, n)
    for i in range(n):
        for j in range(n):
            H[i, j] = G_dense[i][j]
    return H


def eig(H):
    E, Q = mp.eigsy(H)
    n = H.rows
    lam = [E[j] for j in range(n)]
    order = sorted(range(n), key=lambda j: lam[j])
    lam = [lam[j] for j in order]
    Qs = mp.matrix(n, n)
    for a, j in enumerate(order):
        for i in range(n):
            Qs[i, a] = Q[i, j]
    return lam, Qs


def residuals(H, lam, Q):
    n = H.rows
    R = H * Q
    res = mp.mpf(0)
    for j in range(n):
        for i in range(n):
            res = max(res, abs(R[i, j] - lam[j] * Q[i, j]))
    O = Q.T * Q
    orth = max(abs(O[i, j] - (1 if i == j else 0)) for i in range(n) for j in range(n))
    return res, orth


def group(lam, tol):
    """Group sorted eigenvalues into distinct eigenvalues (gap < tol => same)."""
    groups, cur = [], [0]
    for j in range(1, len(lam)):
        if lam[j] - lam[j - 1] < tol:
            cur.append(j)
        else:
            groups.append(cur)
            cur = [j]
    groups.append(cur)
    gaps = [lam[groups[a + 1][0]] - lam[groups[a][-1]] for a in range(len(groups) - 1)]
    within = max([lam[g[-1]] - lam[g[0]] for g in groups] + [mp.mpf(0)])
    return groups, (min(gaps) if gaps else mp.inf), within


def coefficients(Q, groups, s, t):
    """c_lambda = <t|Pi_lambda|s> for each distinct eigenvalue."""
    return [mp.fsum(Q[t, j] * Q[s, j] for j in g) for g in groups]


def S_Fbar(c):
    return mp.fsum(abs(x) for x in c), mp.fsum(x * x for x in c)


def block_masses(Q, lam, intervals, x):
    """||P_i e_x|| for each block interval I_i, and the eigenvalue counts per I_i."""
    out, counts = [], []
    for (lo, hi) in intervals:
        idx = [j for j, l in enumerate(lam) if lo <= l <= hi]
        counts.append(len(idx))
        out.append(mp.sqrt(mp.fsum(Q[x, j] ** 2 for j in idx)))
    return out, counts


def F_sampled(c, glam, max_samples=2_000_000, rel_drop=1e-8):
    """Lower estimate of F_max = sup_tau |sum c e^{-i tau lam}|^2.

    Works with c/S (O(1) weights) in float64. Terms with |c| < rel_drop*S are
    dropped and their total mass subtracted, so the returned value is a rigorous
    lower bound for the sampled sup of the full sum (up to float64 rounding).
    Returns (F_lower, meta)."""
    S = mp.fsum(abs(x) for x in c)
    if S == 0:
        return mp.mpf(0), {}
    w = np.array([float(x / S) for x in c])
    lam = np.array([float(l) for l in glam])
    keep = np.abs(w) >= rel_drop
    dropped = float(np.sum(np.abs(w[~keep])))
    w, lam = w[keep], lam[keep]
    lam = lam - lam.mean()
    spread = lam.max() - lam.min() if len(lam) > 1 else 1.0
    dt = math.pi / (8 * max(spread, 1e-12))
    diffs = np.diff(np.sort(lam))
    diffs = diffs[diffs > 1e-9]
    omega_min = diffs.min() if len(diffs) else 1.0
    T = min(200 * 2 * math.pi / omega_min, max_samples * dt)
    best, best_tau = 0.0, 0.0
    chunk = 200_000
    taus_all = np.arange(0.0, T, dt)
    top = []
    for a in range(0, len(taus_all), chunk):
        taus = taus_all[a:a + chunk]
        amp = np.abs(np.exp(-1j * np.outer(taus, lam)) @ w)
        j = np.argpartition(amp, -5)[-5:] if len(amp) > 5 else np.arange(len(amp))
        top += [(amp[i], taus[i]) for i in j]
    top.sort(reverse=True)
    f = lambda tau: -abs(np.exp(-1j * tau * lam) @ w)
    for amp0, tau0 in top[:20]:
        r = minimize_scalar(f, bounds=(max(0.0, tau0 - dt), tau0 + dt), method="bounded",
                            options=dict(xatol=1e-12))
        if -r.fun > best:
            best, best_tau = -r.fun, r.x
    lower = max(0.0, best - dropped)
    F = (S * mp.mpf(lower)) ** 2
    return F, dict(T=T, dt=dt, n_terms=int(keep.sum()), dropped_rel_mass=dropped,
                   tau_star=best_tau, amp_rel=best)
