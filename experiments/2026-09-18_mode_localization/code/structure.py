"""Stiffness matrix of weakly coupled substructures on graded elastic supports (R-B8).

k substructures, each a chain of 4 unit masses joined by internal springs k_in.
Every mass of substructure i rests on a grounding spring kappa_i = kappa_0 + step*i
(graded foundation).  Consecutive substructures are joined end to end by a weak
coupling spring c.  Unit masses, so the eigenproblem is K phi = omega^2 phi.

K_vv = sum of the stiffnesses attached to mass v (internal + coupling + grounding),
K_uv = -(spring between u and v).  K is supported on the path of 4k vertices; the
coupling springs form a matching, so ||K_out|| = c.  Same dict layout as
families.build, so the certificate code applies unchanged.
"""
from fractions import Fraction as Q


def build_structure(k, k_in=Q(1), c=Q(1, 10), kappa0=Q(1), step=Q(6)):
    edges, h, parts = {}, [], []
    for i in range(k):
        V = list(range(4 * i, 4 * i + 4))
        parts.append(V)
        for a in range(3):
            edges[(V[a], V[a + 1])] = -k_in
        if i + 1 < k:
            edges[(V[3], V[3] + 1)] = -c
    n = 4 * k
    comm = [v // 4 for v in range(n)]
    for v in range(n):
        attached = sum(-x for (u, w), x in edges.items() if v in (u, w))
        h.append(attached + kappa0 + step * comm[v])
    adj = [[] for _ in range(n)]
    for (u, v) in edges:
        adj[u].append(v)
        adj[v].append(u)
    bridges = {e: x for e, x in edges.items() if comm[e[0]] != comm[e[1]]}
    return dict(name="structure", k=k, n=n, edges=edges, h=h, parts=parts, comm=comm, adj=adj,
                s=0, t=n - 1, w=c, bridges=bridges, d=max(len(a) for a in adj),
                mu=max(abs(x) for x in edges.values()),
                physical=dict(k_in=str(k_in), c=str(c), kappa0=str(kappa0), step=str(step)))
