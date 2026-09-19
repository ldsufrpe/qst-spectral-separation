"""Shared primitives, built from scratch from the definitions in main.tex.

NOTHING here is imported from the paper (the paper ships no code: experiments/
is empty).  Graphs, Hamiltonians, distances, spectral projectors and every
derived parameter are assembled directly from Definition 2.1 / Definition 3.1
of the manuscript so that a transcription error in the paper cannot be
reproduced here by construction.
"""
import numpy as np
from collections import deque

# ---------------------------------------------------------------- graph layer

def bfs_all_pairs(n, adj):
    """Unweighted graph distance (Definition 2.1: dist is the distance in G)."""
    D = np.full((n, n), np.inf)
    for src in range(n):
        D[src, src] = 0.0
        q = deque([src])
        while q:
            u = q.popleft()
            for v in adj[u]:
                if D[src, v] == np.inf:
                    D[src, v] = D[src, u] + 1
                    q.append(v)
    return D


def adjacency_from_edges(n, edges):
    adj = [[] for _ in range(n)]
    for (u, v, _w) in edges:
        adj[u].append(v)
        adj[v].append(u)
    return adj


def hamiltonian(n, edges, site):
    """H_{uv}=J_{uv} on E, H_{vv}=site[v]; real symmetric, supported on G."""
    H = np.zeros((n, n))
    for (u, v, w) in edges:
        H[u, v] = w
        H[v, u] = w
    for v in range(n):
        H[v, v] = site[v]
    return H

# ------------------------------------------------------------- spectral layer

def distinct_spectral_projectors(H, tol=1e-9):
    """H = sum_lambda lambda Pi_lambda over DISTINCT eigenvalues (Sec. 2)."""
    w, V = np.linalg.eigh(H)
    groups, cur = [], [0]
    for j in range(1, len(w)):
        if abs(w[j] - w[j - 1]) <= tol * max(1.0, abs(w[j])):
            cur.append(j)
        else:
            groups.append(cur)
            cur = [j]
    groups.append(cur)
    out = []
    for gidx in groups:
        Vg = V[:, gidx]
        out.append((float(np.mean(w[gidx])), Vg @ Vg.T.conj()))
    return out


def S_overlap(H, s, t, tol=1e-9):
    """S(s,t) = sum_lambda |<t|Pi_lambda|s>|   (eq. 2.5)."""
    return sum(abs(P[t, s]) for _, P in distinct_spectral_projectors(H, tol))


def Fbar_exact(H, s, t, tol=1e-9):
    """Fbar = sum_lambda |<t|Pi_lambda|s>|^2  (eq. 2.7)."""
    return sum(abs(P[t, s]) ** 2 for _, P in distinct_spectral_projectors(H, tol))


def amplitude(H, s, t, taus):
    """<t|e^{-i tau H}|s> via the spectral decomposition, from scratch."""
    w, V = np.linalg.eigh(H)
    coef = V[t, :] * V[s, :]                       # real, H real symmetric
    return (coef[None, :] * np.exp(-1j * np.outer(taus, w))).sum(axis=1)

# ------------------------------------------------- (delta,g)-separation layer

def split_in_out(H, comm_of):
    """H = H_in + H_out (Section 3): H_in keeps intra-community entries."""
    n = H.shape[0]
    same = (np.array(comm_of)[:, None] == np.array(comm_of)[None, :])
    Hin = np.where(same, H, 0.0)
    return Hin, H - Hin


def separation_parameters(H, parts, comm_of, d, mu):
    """Every parameter of Definition 3.1 / eq. (3.5)-(3.6), measured on H."""
    Hin, Hout = split_in_out(H, comm_of)
    delta = np.linalg.norm(Hout, 2)
    J, W = [], []
    for Vi in parts:
        sub = H[np.ix_(Vi, Vi)]
        ev = np.linalg.eigvalsh(sub)
        J.append((ev.min(), ev.max()))
        W.append(ev.max() - ev.min())
    g = min(max(J[i][0] - J[j][1], J[j][0] - J[i][1])
            for i in range(len(J)) for j in range(i + 1, len(J)))
    eta = delta / (g - delta)
    ghat = (g - 2 * delta) / (2 * d * mu)
    rho = [1 + 2 * (Wi + 2 * delta) / (np.pi * (g - 2 * delta)) for Wi in W]
    return dict(delta=delta, g=g, J=J, W=W, eta=eta, ghat=ghat, rho=rho,
                Hin=Hin, Hout=Hout)

# ------------------------------------------------------- certificate assembly

def certificate(H, parts, comm_of, d, mu, s, t, kappa, Dmat=None, par=None):
    """Sigma_1, Sigma_2, B_eta, B_E and B, rebuilt from eqs. (6.1)-(6.3),
    Cor. 6.2 and Prop. 3.4 -- no formula copied from the paper's tables."""
    if par is None:
        par = separation_parameters(H, parts, comm_of, d, mu)
    n = H.shape[0]
    if Dmat is None:
        Dmat = bfs_all_pairs(n, adjacency_from_edges(
            n, [(u, v, H[u, v]) for u in range(n) for v in range(u + 1, n)
                if H[u, v] != 0.0]))
    eta, ghat = par['eta'], par['ghat']
    c_kappa = 1.0 / (1.0 - kappa)
    gamma = np.log(1.0 + kappa * ghat)
    A = [c_kappa * par['rho'][i] * np.sqrt(len(parts[i])) / np.sqrt(1 - eta ** 2)
         for i in range(len(parts))]

    def L(i, x):
        b = 1.0 if x in parts[i] else eta
        dist_xVi = min(Dmat[x, u] for u in parts[i])
        return min(b, A[i] * np.exp(-gamma * dist_xVi))

    T = [L(i, s) * L(i, t) for i in range(len(parts))]
    Sigma1 = float(sum(T))
    Sigma2 = float(sum(x * x for x in T))
    b_eta = 2 * eta * np.sqrt(1 - eta ** 2) if eta <= 1 / np.sqrt(2) else 1.0
    B_eta = b_eta if comm_of[s] != comm_of[t] else 1.0
    h = np.diag(H)
    M_H = (h.max() - h.min()) + 2 * d * mu
    B_E = np.sqrt(max(0.0, 1 - (h[s] - h[t]) ** 2 / M_H ** 2))
    B = min(1.0, Sigma1, B_eta, B_E)
    return dict(gamma=gamma, A=A, T=T, Sigma1=Sigma1, Sigma2=Sigma2,
                b_eta=b_eta, B_eta=B_eta, B_E=B_E, B=B,
                F_esp=min(1.0, Sigma1 ** 2), B2=B ** 2, par=par, D=Dmat[s, t])

# --------------------------------------------------------- the three families

def family_cliques(k, m=4, w=0.10, Delta=24.0):
    n = m * k
    parts = [list(range(m * i, m * i + m)) for i in range(k)]
    comm_of = [i // m for i in range(n)]
    edges = []
    for i in range(k):                                   # K_m on each module
        for a in range(m):
            for b in range(a + 1, m):
                edges.append((m * i + a, m * i + b, 1.0))
    for i in range(k - 1):                               # bridge last -> first
        edges.append((m * i + m - 1, m * (i + 1), w))
    site = [Delta * (v // m) for v in range(n)]
    return n, parts, comm_of, edges, site, m, 1.0        # d = m, mu = 1


def family_dimer(k, w=0.05):
    n = 2 * k
    parts = [[2 * i, 2 * i + 1] for i in range(k)]
    comm_of = [i // 2 for i in range(n)]
    edges = [(2 * i, 2 * i + 1, 1.0) for i in range(k)]
    edges += [(2 * i + 1, 2 * i + 2, w) for i in range(k - 1)]
    site = []
    for i in range(k):
        site += [24.0 * i, 24.0 * i + 10.0]
    return n, parts, comm_of, edges, site, 2, 1.0


def family_p4(k, w=0.30):
    n = 4 * k
    parts = [list(range(4 * i, 4 * i + 4)) for i in range(k)]
    comm_of = [i // 4 for i in range(n)]
    edges = []
    for i in range(k):
        for a in range(3):
            edges.append((4 * i + a, 4 * i + a + 1, 1.0))
    for i in range(k - 1):
        edges.append((4 * i + 3, 4 * (i + 1), w))
    site = [5.0 * (v // 4) for v in range(n)]
    return n, parts, comm_of, edges, site, 2, 1.0


def assemble(family):
    n, parts, comm_of, edges, site, d, mu = family
    H = hamiltonian(n, edges, site)
    Dmat = bfs_all_pairs(n, adjacency_from_edges(n, edges))
    return n, parts, comm_of, H, Dmat, d, mu
