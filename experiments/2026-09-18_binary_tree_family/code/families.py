"""The three example families of Section 7 (sec:exemplos), built from the text.

Every entry is an exact integer or rational. Nothing here is floating point, so the
same object feeds the Arb certification and, as float64 or mpmath, the
diagonalization benchmarks of the sibling experiments.
"""
from collections import deque
from fractions import Fraction as Q


def _family(name, k):
    """Return (n, edges, h, parts, s, t, w) for family `name` with k modules.

    edges: dict {(u, v): weight} with u < v, internal and bridge edges together.
    h: list of site energies (exact). parts: list of vertex lists V_0..V_{k-1}.
    """
    edges, h, parts = {}, [], []
    if name == "cliques":          # Example ex:cliques: K4 modules, w=0.10, detuning 24
        m, w = 4, Q(1, 10)
        for i in range(k):
            V = list(range(m * i, m * i + m))
            parts.append(V)
            for a in range(m):
                for b in range(a + 1, m):
                    edges[(V[a], V[b])] = Q(1)
            h += [Q(24 * i)] * m
            if i + 1 < k:
                edges[(V[-1], V[-1] + 1)] = w
        s, t = 0, m * k - 1
    elif name == "dimer":          # Example ex:dimero: H_i = [[24i,1],[1,24i+10]], w=0.05
        w = Q(1, 20)
        for i in range(k):
            V = [2 * i, 2 * i + 1]
            parts.append(V)
            edges[(V[0], V[1])] = Q(1)
            h += [Q(24 * i), Q(24 * i + 10)]
            if i + 1 < k:
                edges[(V[1], V[1] + 1)] = w
        s, t = 0, 2 * k - 1
    elif name == "p4":             # Example ex:p4: P4 modules, w=0.3, site energy 5i
        m, w = 4, Q(3, 10)
        for i in range(k):
            V = list(range(m * i, m * i + m))
            parts.append(V)
            for a in range(m - 1):
                edges[(V[a], V[a + 1])] = Q(1)
            h += [Q(5 * i)] * m
            if i + 1 < k:
                edges[(V[-1], V[-1] + 1)] = w
        s, t = 0, m * k - 1
    else:
        raise ValueError(name)
    n = len(h)
    return n, edges, h, parts, s, t, w


def build(name, k):
    n, edges, h, parts, s, t, w = _family(name, k)
    comm = [None] * n
    for i, V in enumerate(parts):
        for v in V:
            comm[v] = i
    adj = [[] for _ in range(n)]
    for (u, v) in edges:
        adj[u].append(v)
        adj[v].append(u)
    bridges = {e: x for e, x in edges.items() if comm[e[0]] != comm[e[1]]}
    return dict(name=name, k=k, n=n, edges=edges, h=h, parts=parts, comm=comm,
                adj=adj, s=s, t=t, w=w, bridges=bridges,
                d=max(len(a) for a in adj),
                mu=max(abs(x) for x in edges.values()))


def bfs(adj, src):
    dist = [-1] * len(adj)
    dist[src] = 0
    q = deque([src])
    while q:
        u = q.popleft()
        for v in adj[u]:
            if dist[v] < 0:
                dist[v] = dist[u] + 1
                q.append(v)
    return dist


def dist_to_parts(G, x):
    dx = bfs(G["adj"], x)
    return [min(dx[v] for v in V) for V in G["parts"]], dx


def bridges_form_matching(G):
    seen = set()
    for (u, v) in G["bridges"]:
        if u in seen or v in seen:
            return False
        seen.update((u, v))
    return True


def block_matrix(G, i):
    """Integer/rational block H_i = H|_{V_i} as a list of lists (exact)."""
    V = G["parts"][i]
    idx = {v: a for a, v in enumerate(V)}
    M = [[Q(0)] * len(V) for _ in V]
    for v in V:
        M[idx[v]][idx[v]] = G["h"][v]
    for (u, v), x in G["edges"].items():
        if u in idx and v in idx:
            M[idx[u]][idx[v]] = M[idx[v]][idx[u]] = x
    return M


def dense(G, conv=float):
    n = G["n"]
    H = [[conv(0)] * n for _ in range(n)]
    for v in range(n):
        H[v][v] = conv(G["h"][v])
    for (u, v), x in G["edges"].items():
        H[u][v] = H[v][u] = conv(x)
    return H
