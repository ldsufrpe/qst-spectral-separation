"""Binary tree of K4 modules (finding R-B4): a non-chain family for Theorem C / (GR).

Module labels l = 0..k-1 in BFS order (children of l: 2l+1, 2l+2), k = 2^{r+1}-1.
Module l occupies vertices 4l..4l+3 and induces K4 (coupling 1).
Bridges (weight w): parent vertex 4p+1 -> first child's vertex 4c, parent vertex
4p+2 -> second child's vertex 4c.  Every vertex carries at most one bridge, so the
bridges form a matching and ||H_out|| = w.  Site energy 24*l on module l, exactly
the blocks of Example 7.1, so W = 4, g = 20.
s = vertex 3 of the root (no bridge); t = vertex 3 of the leftmost leaf.
"""
from fractions import Fraction as Q


def build_tree(r, w=Q(1, 10), Delta=Q(24)):
    k = 2 ** (r + 1) - 1
    edges, h, parts = {}, [], []
    for l in range(k):
        V = list(range(4 * l, 4 * l + 4))
        parts.append(V)
        for a in range(4):
            for b in range(a + 1, 4):
                edges[(V[a], V[b])] = Q(1)
        h += [Delta * l] * 4
        for slot, c in ((1, 2 * l + 1), (2, 2 * l + 2)):
            if c < k:
                edges[(4 * l + slot, 4 * c)] = w
    n = 4 * k
    comm = [v // 4 for v in range(n)]
    adj = [[] for _ in range(n)]
    for (u, v) in edges:
        adj[u].append(v)
        adj[v].append(u)
    bridges = {e: x for e, x in edges.items() if comm[e[0]] != comm[e[1]]}
    leftmost_leaf = 2 ** r - 1
    return dict(name=f"tree_r{r}", r=r, k=k, n=n, edges=edges, h=h, parts=parts, comm=comm,
                adj=adj, s=3, t=4 * leftmost_leaf + 3, w=w, bridges=bridges,
                d=max(len(a) for a in adj), mu=max(abs(x) for x in edges.values()))
