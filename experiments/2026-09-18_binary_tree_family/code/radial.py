"""Exact radial counts for (GR) and C_theta (Theorem cor:radial, eq:Ctheta, eq:GR).

For every vertex x, r_i(x) = dist(x, V_i).  N_x(R) = #{i : r_i(x) <= R} is a right-
continuous step function and e^{-aR} is decreasing, so
    sup_R N_x(R) e^{-aR} = max over the jump points R in {r_i(x)}.
Hence C(a) := max_x sup_R N_x(R) e^{-aR} is computed exactly (no sampling in R).
C_theta := max_x sum_i e^{-theta r_i(x)} is computed exactly as well.
All-pairs BFS is done with scipy.sparse.csgraph in source batches.
"""
import math

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import shortest_path


def module_distance_profiles(G, batch=512):
    """Yield, for batches of vertices x, the int array r[x, i] = dist(x, V_i)."""
    n, k = G["n"], G["k"]
    rows, cols = [], []
    for (u, v) in G["edges"]:
        rows += [u, v]
        cols += [v, u]
    A = csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(n, n))
    comm = np.array(G["comm"])
    for start in range(0, n, batch):
        src = np.arange(start, min(n, start + batch))
        D = shortest_path(A, unweighted=True, indices=src, directed=False)
        Di = D.astype(np.int64)
        m = n // k
        if all(G["parts"][i] == list(range(m * i, m * i + m)) for i in range(k)):
            R = Di.reshape(len(src), k, m).min(axis=2)          # contiguous equal modules
        else:
            R = np.stack([Di[:, comm == i].min(axis=1) for i in range(k)], axis=1)
        yield src, R


def radial_constants(G, a_grid, theta_grid):
    """Return dict with C(a) for a in a_grid and C_theta for theta in theta_grid,
    plus the argmax vertex for each."""
    Ca = {a: (0.0, None) for a in a_grid}
    Ct = {th: (0.0, None) for th in theta_grid}
    for src, R in module_distance_profiles(G):
        Rs = np.sort(R, axis=1)                       # r_(1) <= r_(2) <= ...
        counts = np.arange(1, R.shape[1] + 1)         # N_x at the j-th jump (ties: last wins)
        for a in a_grid:
            val = (counts[None, :] * np.exp(-a * Rs)).max(axis=1)
            j = int(val.argmax())
            if val[j] > Ca[a][0]:
                Ca[a] = (float(val[j]), int(src[j]))
        for th in theta_grid:
            val = np.exp(-th * R).sum(axis=1)
            j = int(val.argmax())
            if val[j] > Ct[th][0]:
                Ct[th] = (float(val[j]), int(src[j]))
    return Ca, Ct
