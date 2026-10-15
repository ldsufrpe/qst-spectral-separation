"""Rigorous (Arb) Neumann-route certificates on the block masses ||P_i e_x||.

Formulas: derivacoes/2026-10-08_neumann_pontes/derivacao.md, sections 2-6
(Lemma N0, Theorem N1, Proposition N2, Corollary N2', Proposition N3, diagnostic N3x).

Ported from experiments/toy/2026-10-08_neumann_calibration/code/run.py (SHA-256 logged
in experiment_log.json): rig_setup, bounds_123, stadium_segments and the dense N3x
(here `n3x_upper_dense`, used only by crosscheck.py). The toy's dense arb_mat.solve and
semipositivity check are replaced by mlu_solve (componentwise-stable M-matrix LU).
New here: mlu_solve, n2prime (three forms), ladder_certified, community_path,
n3x_upper_reduced, and resolvent-identity segment enclosures (_U_block).

Every parameter that enters a bound is an exact rational on the safe side
(g, g_ij lower; W, eta upper; delta exact). Every per-block bound is returned as an
Arb ball whose upper endpoint is a rigorous upper bound on ||P_i e_x||.
"""
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
for p in (ROOT / "experiments/2026-09-18_exact_benchmarks/code",
          ROOT / "experiments/2026-09-18_interval_certification/code"):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

from flint import arb, acb, acb_mat, fmpq                 # noqa: E402

from families import build, bfs                           # noqa: E402
from certify import block_spectrum, upper_q, lower_q, amin    # noqa: E402

KAPPAS = (fmpq(1, 2), fmpq(7, 10), fmpq(4, 5), fmpq(9, 10), fmpq(19, 20))


# ------------------------------------------------------------------ utilities
def fq(x):
    """Fraction -> fmpq (exact)."""
    return fmpq(x.numerator, x.denominator)


def lg10(x):
    """log10 of the upper endpoint of x (an arb or fmpq); -inf for 0."""
    u = arb(x).upper()
    return float(u.log() / arb(10).log()) if lower_q(u) > 0 else -math.inf


class NotMMatrix(Exception):
    """I - A could not be certified a nonsingular M-matrix (a pivot not provably > 0)."""


def mlu_solve(n, A, cols):
    """Solve (I - A) Y = C for a non-negative sparse A ({row: {col: arb >= 0}}), by
    Gaussian elimination without pivoting in index order. Returns Y as a list of rows.

    Certificate: I - A is a Z-matrix; all pivots certified > 0 means all leading principal
    minors are positive, so I - A is a nonsingular M-matrix (Berman-Plemmons, Thm 6.2.3),
    i.e. rho(A) < 1 and (I - A)^{-1} >= 0. Raises NotMMatrix otherwise.
    Componentwise stable in ball arithmetic: off-diagonal updates and both substitutions
    only add non-negative terms, so components of size 1e-270 keep their relative accuracy
    (a dense arb_mat.solve bounds the error normwise and would swamp them)."""
    zero = arb(0)
    rows = [{i: arb(1)} for i in range(n)]
    colnz = [{i} for i in range(n)]
    for i, d in A.items():
        for j, a in d.items():
            rows[i][j] = rows[i].get(j, zero) - a
            colnz[j].add(i)
    m = len(cols)
    b = [[cols[c][i] for c in range(m)] for i in range(n)]
    for kk in range(n):
        piv = rows[kk][kk]
        if not lower_q(piv) > 0:
            raise NotMMatrix(kk)
        upper_kk = [(j, a) for j, a in rows[kk].items() if j > kk]
        for i in sorted(r for r in colnz[kk] if r > kk):
            lik = rows[i].pop(kk, None)
            if lik is None:
                continue
            lik = lik / piv
            for j, akj in upper_kk:
                rows[i][j] = rows[i].get(j, zero) - lik * akj
                colnz[j].add(i)
            bk = b[kk]
            bi = b[i]
            for c in range(m):
                bi[c] = bi[c] - lik * bk[c]
    y = [None] * n
    for kk in range(n - 1, -1, -1):
        acc = list(b[kk])
        for j, akj in rows[kk].items():
            if j > kk:
                for c in range(m):
                    acc[c] = acc[c] - akj * y[j][c]
        piv = rows[kk][kk]
        y[kk] = [v / piv for v in acc]
    return y


def _unit(n, idx):
    return [arb(1) if u == idx else arb(0) for u in range(n)]


def _bridge_list(P):
    """Directed bridge entries (u, u', |w|) of |V|."""
    return [(u, v, P["absV"][u][v]) for u in range(P["n"]) for v in range(P["n"])
            if P["absV"][u][v] != 0]


def catalan_B(z):
    return (1 - (1 - 4 * z).sqrt()) / (2 * z)


def dfact(m):
    return math.prod(range(m, 0, -2)) if m > 0 else 1


# ------------------------------------------------------------ rigorous set-up
def rig_setup(name, k):
    """Exact/rigorous data of one instance of the three Section 7 families."""
    G = build(name, k)
    n, parts, comm = G["n"], G["parts"], G["comm"]
    lo, hi = [], []
    for i in range(k):
        _, _, roots = block_spectrum(G, i)
        vals = [r for r, _ in roots]
        mn = mx = vals[0]
        for v in vals[1:]:
            mn = amin(mn, v, "lo")
            mx = -amin(-mx, -v, "hi")
        lo.append(mn)
        hi.append(mx)
    lo_q = [lower_q(x) for x in lo]          # rational outer enclosure of J_i
    hi_q = [upper_q(x) for x in hi]
    W = [hi_q[i] - lo_q[i] for i in range(k)]
    gij = [[None] * k for _ in range(k)]
    for i in range(k):
        for j in range(k):
            if i != j:
                gij[i][j] = max(lo_q[j] - hi_q[i], lo_q[i] - hi_q[j])
    g = min(gij[i][j] for i in range(k) for j in range(k) if i != j)
    seen = set()
    for (u, v) in G["bridges"]:
        assert u not in seen and v not in seen, "bridges are not a matching"
        seen.update((u, v))
    delta = max(abs(fq(x)) for x in G["bridges"].values())   # = ||H_out|| (matching)
    assert g > 2 * delta, "(H2) fails rigorously"
    rho = g / 2
    eta = delta / (g - delta)
    absV = [[fmpq(0)] * n for _ in range(n)]
    Vblk = [[fmpq(0)] * k for _ in range(k)]
    for (u, v), x in G["bridges"].items():
        absV[u][v] = absV[v][u] = abs(fq(x))
        a, b = comm[u], comm[v]
        # single bridge edge per community pair in these families => ||Q_a V Q_b|| = |w|
        Vblk[a][b] = Vblk[b][a] = max(Vblk[a][b], abs(fq(x)))
    npairs = {}
    for (u, v) in G["bridges"]:
        key = tuple(sorted((comm[u], comm[v])))
        npairs[key] = npairs.get(key, 0) + 1
    assert all(c == 1 for c in npairs.values()), "multi-edge community pair: Vblk needs a norm"
    cadj = [[l for l in range(k) if Vblk[j][l] > 0] for j in range(k)]
    Ncomm = [bfs(cadj, a) for a in range(k)]
    inner = []
    for V in parts:
        idx = {v: a for a, v in enumerate(V)}
        adj = [[idx[u] for u in G["adj"][v] if u in idx] for v in V]
        dmat = [bfs(adj, a) for a in range(len(V))]
        dj = max(len(a) for a in adj)
        muj = max([abs(fq(x)) for (u, v), x in G["edges"].items()
                   if u in idx and v in idx] + [fmpq(0)])
        inner.append((dmat, dj, muj))
    Hblocks = []
    for V in parts:
        idx = {v: a for a, v in enumerate(V)}
        Hj = [[fmpq(0)] * len(V) for _ in V]
        for v in V:
            Hj[idx[v]][idx[v]] = fq(G["h"][v])
        for (u, v), x in G["edges"].items():
            if u in idx and v in idx:
                Hj[idx[u]][idx[v]] = Hj[idx[v]][idx[u]] = fq(x)
        Hblocks.append(Hj)
    B = sorted(seen)                                       # bridge endpoints ("doors")
    D = bfs(G["adj"], G["s"])[G["t"]]
    return dict(G=G, name=name, k=k, n=n, parts=parts, comm=comm, lo=lo_q, hi=hi_q, W=W,
                gij=gij, g=g, delta=delta, rho=rho, eta=eta, absV=absV, Vblk=Vblk,
                cadj=cadj, Ncomm=Ncomm, inner=inner, Hblocks=Hblocks, B=B, D=D,
                s=G["s"], t=G["t"])


def community_path(P):
    """The community graph is a path (Corollary N2' hypothesis): connected, k-1 edges,
    degrees <= 2. Returns the path order of the communities, or None."""
    k, cadj = P["k"], P["cadj"]
    if sum(len(a) for a in cadj) != 2 * (k - 1) or any(len(a) > 2 for a in cadj):
        return None
    if any(d < 0 for d in P["Ncomm"][0]):
        return None
    ends = [j for j in range(k) if len(cadj[j]) <= 1]
    order = sorted(range(k), key=lambda j: P["Ncomm"][ends[0]][j])
    return order


def ladder_certified(P):
    """g_il >= N(i,l) g for all pairs (exact rationals): the ordered-ladder condition."""
    k = P["k"]
    return all(P["gij"][i][l] >= P["Ncomm"][i][l] * P["g"]
               for i in range(k) for l in range(k) if i != l)


# ------------------------------------------------------------ N1, N2, N3
def bounds_123(P, kappas=KAPPAS):
    """Per-block bounds {cert: {(i, x): arb}} for N1, N2, N3 and N3vN2 = min(N3, N2),
    plus certification flags (rho < 1 certified for every walk matrix used)."""
    k, n = P["k"], P["n"]
    rho, delta, eta = arb(P["rho"]), arb(P["delta"]), arb(P["eta"])
    q = delta / rho
    sq = (1 - eta ** 2).sqrt()
    xs = (P["s"], P["t"])
    a_of = {x: P["comm"][x] for x in xs}
    out = {key: {} for key in ("N1", "N2", "N3", "N3vN2")}
    cert = {"N2": True, "N3": True}
    bridges = _bridge_list(P)
    blk_of = P["comm"]
    for i in range(k):
        pref = (2 * arb(P["W"][i]) + 2 * arb.pi() * rho) / (2 * arb.pi())
        dhat = [1 / rho if j == i else 1 / (arb(P["gij"][i][j]) - rho) for j in range(k)]
        T = {j: {l: arb(P["Vblk"][j][l]) * dhat[l] for l in range(k) if P["Vblk"][j][l] != 0}
             for j in range(k)}
        try:
            sol2 = mlu_solve(k, T, [_unit(k, a_of[x]) for x in xs])
            ok2 = True
        except NotMMatrix:
            sol2, ok2 = None, False
        cert["N2"] &= ok2
        best3 = {x: None for x in xs}
        for kap in kappas:
            ck = 1 / (1 - arb(kap))
            Rh = {}                                   # (u, v) -> arb, u, v in the same block
            for j, V in enumerate(P["parts"]):
                dj = rho if j == i else arb(P["gij"][i][j]) - rho
                dmat, degj, muj = P["inner"][j]
                for aa, u in enumerate(V):
                    for bb, v in enumerate(V):
                        dd = dmat[aa][bb]
                        if dd < 0:
                            continue
                        if degj * muj > 0:
                            ct = ck / dj * (1 + arb(kap) * dj / (degj * arb(muj))) ** (-dd)
                            Rh[(u, v)] = amin(1 / dj, ct, "Rh")
                        elif dd == 0:
                            Rh[(u, v)] = 1 / dj
            T3 = {}
            for (u, up, w) in bridges:                # (|V| Rh)_{uv} = sum_u' |w| Rh_{u'v}
                row = T3.setdefault(u, {})
                for v in P["parts"][blk_of[up]]:
                    if (up, v) in Rh:
                        row[v] = row.get(v, arb(0)) + arb(w) * Rh[(up, v)]
            try:
                y = mlu_solve(n, T3, [_unit(n, x) for x in xs])
            except NotMMatrix:
                cert["N3"] = False
                continue
            for c, x in enumerate(xs):
                if a_of[x] == i:
                    continue
                ssq = arb(0)
                for u in P["parts"][i]:
                    mu = sum((Rh[(u, v)] * y[v][c] for v in P["parts"][i] if (u, v) in Rh), arb(0))
                    ssq += mu ** 2
                val = ssq.sqrt()
                best3[x] = val if best3[x] is None else amin(best3[x], val, "kappa")
        for c, x in enumerate(xs):
            a = a_of[x]
            if a == i:
                for key in out:
                    out[key][(i, x)] = arb(1)
                continue
            N = P["Ncomm"][a][i]
            geo = q ** N / (1 - q) / rho
            out["N1"][(i, x)] = amin(eta, pref * geo / sq, "N1")
            n2 = pref * amin(dhat[i] * sol2[i][c], geo, "N2") / sq if ok2 else pref * geo / sq
            out["N2"][(i, x)] = amin(eta, n2, "N2eta")
            if best3[x] is not None:
                n3 = amin(eta, pref * amin(best3[x], geo, "N3") / sq, "N3eta")
            else:
                n3 = out["N1"][(i, x)]
                cert["N3"] = False
            out["N3"][(i, x)] = n3
            out["N3vN2"][(i, x)] = amin(n3, out["N2"][(i, x)], "N3vN2")
    return out, cert


# ------------------------------------------------------------ N2' (chain corollary)
def n2prime(P):
    """Corollary N2' on a chain of communities, three forms (derivacao.md sec. 4):
      stated : global delta, product over the sorted distances G_(m) of block i;
      path   : the proof's own first-arrival product over the visited blocks, with the
               coupling of each bridge crossed (valid on any chain; <= stated);
      ladder : global delta, q^N/(2N-1)!!, only where the ladder condition is certified.
    Returns ({form: {(i, x): arb}}, info). Forms whose walk-count series diverges
    (4 z >= 1) are omitted."""
    order = community_path(P)
    assert order is not None, "community graph is not a path"
    k = P["k"]
    rho, delta, eta = arb(P["rho"]), arb(P["delta"]), arb(P["eta"])
    sq = (1 - eta ** 2).sqrt()
    q = delta / rho
    dmax = arb(max(P["Vblk"][j][l] for j in range(k) for l in range(k)))
    qb = dmax / rho
    lad = ladder_certified(P)
    xs = (P["s"], P["t"])
    forms = {"N2p_stated": {}, "N2p_path": {}}
    if lad:
        forms["N2p_ladder"] = {}
    ok = {"N2p_stated": lower_q(1 - 4 * q ** 2) > 0, "N2p_path": lower_q(1 - 4 * qb ** 2) > 0}
    ok["N2p_ladder"] = ok["N2p_stated"]
    pos = {c: r for r, c in enumerate(order)}
    for i in range(k):
        rt = 1 + 2 * arb(P["W"][i]) / (arb.pi() * 2 * rho)     # rho~_i = 1 + 2W/(pi g)
        x_of = {j: 1 / (arb(P["gij"][i][j]) - rho) for j in range(k) if j != i}
        G_sorted = sorted((P["gij"][i][l] for l in range(k) if l != i))
        for x in xs:
            a = P["comm"][x]
            if a == i:
                for f in forms:
                    forms[f][(i, x)] = arb(1)
                continue
            N = P["Ncomm"][a][i]
            # stated
            if ok["N2p_stated"]:
                prod = arb(1)
                for m in range(N):
                    prod *= delta / (arb(G_sorted[m]) - rho)
                val = rt / sq * catalan_B(q ** 2) ** N / (1 - 4 * q ** 2).sqrt() * prod
                forms["N2p_stated"][(i, x)] = amin(eta, val, "N2p_stated")
            # path: blocks a = j0, j1, ..., jN = i along the chain
            if ok["N2p_path"]:
                step = 1 if pos[i] > pos[a] else -1
                seq = [order[r] for r in range(pos[a], pos[i] + step, step)]
                prod = arb(1)
                for l, lnext in zip(seq[:-1], seq[1:]):
                    prod *= x_of[l] * arb(P["Vblk"][l][lnext])
                val = rt / sq * catalan_B(qb ** 2) ** N / (1 - 4 * qb ** 2).sqrt() * prod
                forms["N2p_path"][(i, x)] = amin(eta, val, "N2p_path")
            if lad and ok["N2p_ladder"]:
                val = rt / sq * catalan_B(q ** 2) ** N / (1 - 4 * q ** 2).sqrt() \
                    * q ** N / dfact(2 * N - 1)
                forms["N2p_ladder"][(i, x)] = amin(eta, val, "N2p_ladder")
    forms = {f: v for f, v in forms.items() if ok.get(f, False) and v}
    return forms, dict(ladder_certified=lad, q=str(q), q_bridge=str(qb))


# ------------------------------------------------------------ N3x (rigorous)
def stadium_segments(lo, hi, rho, M):
    """Segments of the boundary of the rho-neighbourhood of [lo, hi] (all fmpq).
    Returns ([(z0, r, ln)], per): z0 a point of the curve (tiny acb ball), r an upper
    bound on max |z - z0| over the segment, ln an upper bound on its length."""
    L1 = arb(hi - lo)
    per = 2 * L1 + 2 * arb.pi() * arb(rho)
    ns = max(1, int(round(M * float(L1 / per))))
    na = max(1, int(round(M * float(arb.pi() * arb(rho) / per))))
    segs = []
    for sign in (-1, 1):
        for a in range(ns):
            x0 = lo + (hi - lo) * fmpq(a, ns)
            x1 = lo + (hi - lo) * fmpq(a + 1, ns)
            segs.append((acb(arb((x0 + x1) / 2), arb(sign * rho)), arb((x1 - x0) / 2),
                         arb(x1 - x0)))
    dth = arb.pi() / na
    for c, th0 in ((hi, -arb.pi() / 2), (lo, arb.pi() / 2)):
        for a in range(na):
            th = th0 + (a + arb(1) / 2) * dth
            z0 = acb(arb(c) + arb(rho) * th.cos(), arb(rho) * th.sin())
            segs.append((z0, arb(rho) * dth / 2, arb(rho) * dth))   # arc half-length >= |z-z0|
    return segs, per


def _dmin(P, i):
    """Lower bounds on dist(z, J_j) for z on Gamma_i (exact rationals): rho for j = i,
    gap_ij - rho otherwise (gap between the rational outer enclosures)."""
    out = []
    for j in range(P["k"]):
        if j == i:
            out.append(P["rho"])
        else:
            out.append(max(P["lo"][j] - P["hi"][i], P["lo"][i] - P["hi"][j]) - P["rho"])
    return out


def _U_block(Hb, Ib, z0, r, dmin, K=4):
    """Entrywise upper bounds of |(H_j - z)^{-1}| over a segment {|z - z0| <= r}.

    Taylor expansion at z0: R(z) = sum_m (z - z0)^m R(z0)^{m+1}, convergent for
    r ||R(z0)|| < 1. H_j is normal and dist(z0, spec H_j) >= dmin, so ||R(z0)|| <= 1/dmin and
        |R(z)_uv| <= sum_{m<K} r^m |(R(z0)^{m+1})_uv| + (r/dmin)^K / dmin / (1 - r/dmin).
    The first K terms are entrywise, so entries that are small by cancellation inside the
    block (the door-to-door "door factor", ~1/d^4 across a P4 module) stay small. An
    absolute pad r/dmin^2 (first version, pass 1) swamped them on the P4 rows."""
    R0 = (Hb - Ib * z0).inv()
    n = Hb.nrows()
    ru = r.upper()
    d = arb(dmin)
    assert nm_upper(ru / d) < 1, "segment radius not below the distance to the spectrum"
    acc = [[abs(R0[a, b]) for b in range(n)] for a in range(n)]
    Pw = R0
    for m in range(1, K):
        Pw = Pw * R0
        rm = ru ** m
        for a in range(n):
            for b in range(n):
                acc[a][b] += rm * abs(Pw[a, b])
    tail = (ru / d) ** K / d / (1 - ru / d)
    return [[(acc[a][b] + tail).upper() for b in range(n)] for a in range(n)]


def nm_upper(x):
    return upper_q(arb(x))


def _riesz_checks(P, i):
    lo, hi, rho, delta = P["lo"][i], P["hi"][i], P["rho"], P["delta"]
    assert rho > delta, "contour does not enclose I_i"
    for j in range(P["k"]):
        if j != i:
            gap = max(P["lo"][j] - hi, lo - P["hi"][j])
            assert gap - rho > delta, "contour too close to another block"


def _block_mats(P):
    Hb = [acb_mat([[acb(arb(x)) for x in row] for row in Hj]) for Hj in P["Hblocks"]]
    Ib = [acb_mat([[acb(1) if a == b else acb(0) for b in range(len(Hj))]
                   for a in range(len(Hj))]) for Hj in P["Hblocks"]]
    return Hb, Ib


def _finish(P, i, acc):
    eta = arb(P["eta"])
    sq = (1 - eta ** 2).sqrt()
    res = {}
    for c, x in enumerate((P["s"], P["t"])):
        res[x] = arb(1) if P["comm"][x] == i else amin(eta, acc[c] / (2 * arb.pi()) / sq, "N3x")
    return res


def _n3x_U(P, Hb, Ib, z0, r, dmin):
    """Per-block entrywise upper bounds of |R_j(z)| over the segment, as {(u, v): arb}."""
    U = {}
    for j, V in enumerate(P["parts"]):
        Uj = _U_block(Hb[j], Ib[j], z0, r, dmin[j])
        for aa, u in enumerate(V):
            for bb, v in enumerate(V):
                U[(u, v)] = Uj[aa][bb]
    return U


def n3x_upper_dense(P, i, M):
    """Rigorous N3x upper sum on the full n x n system (crosscheck of the reduction)."""
    _riesz_checks(P, i)
    n = P["n"]
    segs, per = stadium_segments(P["lo"][i], P["hi"][i], P["rho"], M)
    xs = (P["s"], P["t"])
    Hb, Ib = _block_mats(P)
    bridges = _bridge_list(P)
    dmin = _dmin(P, i)
    acc = [arb(0), arb(0)]
    for z0, r, ln in segs:
        U = _n3x_U(P, Hb, Ib, z0, r, dmin)
        T = {}
        for (u, up, w) in bridges:                    # sum over all bridge partners u'
            row = T.setdefault(u, {})
            for v in P["parts"][P["comm"][up]]:
                row[v] = row.get(v, arb(0)) + arb(w) * U[(up, v)]
        try:
            y = mlu_solve(n, T, [_unit(n, x) for x in xs])
        except NotMMatrix:
            raise RuntimeError(f"rho(|V|U) < 1 not certified (block {i}, M={M})")
        for c, x in enumerate(xs):
            if P["comm"][x] == i:
                continue
            ssq = arb(0)
            for u in P["parts"][i]:
                ssq += sum((U[(u, v)] * y[v][c] for v in P["parts"][i]), arb(0)) ** 2
            acc[c] += ln.upper() * ssq.sqrt().upper()
    return _finish(P, i, acc), len(segs)


def n3x_upper_reduced(P, i, M):
    """Rigorous N3x upper sum, reduced to the bridge endpoints B.

    With w = U y and y = e_x + |V| w, only y_B is unknown:
        (I - |V|_BB U_BB) y_B = (e_x)_B + [x not in B] |V|_BB U_Bx,
        (Q_i m)_u = sum_{v in B cap V_i} U_uv y_v      (i != comm(x)).
    spec(|V| U) = spec(|V|_BB U_BB) + zeros, so the M-matrix certificate of the reduced
    system (mlu_solve) gives rho(|V|U) < 1. U holds entrywise upper bounds of |R_j(z)|
    over each whole segment (_U_block); by monotonicity of the non-negative Neumann
    series the result bounds the integrand at every z of the segment."""
    _riesz_checks(P, i)
    B = P["B"]
    nb = len(B)
    pos = {u: a for a, u in enumerate(B)}
    segs, per = stadium_segments(P["lo"][i], P["hi"][i], P["rho"], M)
    xs = (P["s"], P["t"])
    Hb, Ib = _block_mats(P)
    bridges = _bridge_list(P)
    doors = [[v for v in V if v in pos] for V in P["parts"]]
    dmin = _dmin(P, i)
    acc = [arb(0), arb(0)]
    zero = arb(0)
    for z0, r, ln in segs:
        U = _n3x_U(P, Hb, Ib, z0, r, dmin)
        A = {}
        cols = [[zero] * nb for _ in xs]
        for (u, up, w) in bridges:                    # sum over all bridge partners u'
            ju = P["comm"][up]
            row = A.setdefault(pos[u], {})
            for v in doors[ju]:
                row[pos[v]] = row.get(pos[v], zero) + arb(w) * U[(up, v)]
            for c, x in enumerate(xs):
                if x not in pos and P["comm"][x] == ju:
                    cols[c][pos[u]] = cols[c][pos[u]] + arb(w) * U[(up, x)]
        for c, x in enumerate(xs):
            if x in pos:
                cols[c][pos[x]] = arb(1)
        try:
            y = mlu_solve(nb, A, cols)
        except NotMMatrix:
            raise RuntimeError(f"rho(|V|U) < 1 not certified (block {i}, M={M})")
        for c, x in enumerate(xs):
            if P["comm"][x] == i:
                continue
            ssq = zero
            for u in P["parts"][i]:
                ssq += sum((U[(u, v)] * y[pos[v]][c] for v in doors[i]), zero) ** 2
            acc[c] += ln.upper() * ssq.sqrt().upper()
    return _finish(P, i, acc), len(segs)
