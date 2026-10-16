"""R-C3: the chain length k at which each certificate guarantees Sigma_1 < sqrt2 - 1, so that
F_Bloch < 2/3 at every time (cor:fent), next to the exact overlap sum S (see ../plan.md).

Read-only imports (the files are never edited):
  experiments/2026-10-08_neumann_14rows/code/neumann.py
      rig_setup, bounds_123 (N1 = bridge, N2 = walk-sum), n2prime (N2p_path = chain),
      n3x_upper_reduced (contour), fq, lg10
  experiments/2026-09-18_interval_certification/code/certify.py
      parameters, geometry, sigma1 (Combes-Thomas baseline), upper_q, lower_q, amin
  experiments/2026-09-18_exact_benchmarks/code/families.py      build, dense
  experiments/2026-09-18_exact_benchmarks/code/spectral_mp.py   to_mp, eig, residuals, group,
      coefficients, S_Fbar, block_masses
Copied, not imported: mpq and the adaptive-precision rule of
experiments/2026-09-18_exact_benchmarks/code/run_benchmark.py (l.31, l.102-124); importing that
module would pull in its own logutil under the same module name.

neumann.py resolves `families` to the copy in interval_certification, which is identical to the
one in exact_benchmarks (both SHA-256 are logged). Stronger bridges enter through build_w, which
rig_setup sees through a scoped rebinding of the module-level name neumann.build.
"""
import contextlib
import math
import sys
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
NEU = ROOT / "experiments/2026-10-08_neumann_14rows/code"
if str(NEU) not in sys.path:
    sys.path.insert(0, str(NEU))

import mpmath as mp                                     # noqa: E402
from flint import arb, fmpq                             # noqa: E402

import neumann as nm                                    # noqa: E402  (inserts the two paths below)
import certify as C                                     # noqa: E402
import families as FAM                                  # noqa: E402
import spectral_mp as SP                                # noqa: E402

REUSED = [NEU / "neumann.py",
          ROOT / "experiments/2026-09-18_interval_certification/code/certify.py",
          ROOT / "experiments/2026-09-18_interval_certification/code/families.py",
          ROOT / "experiments/2026-09-18_exact_benchmarks/code/families.py",
          ROOT / "experiments/2026-09-18_exact_benchmarks/code/spectral_mp.py",
          ROOT / "experiments/2026-09-18_exact_benchmarks/code/run_benchmark.py",
          NEU / "logutil.py"]

# ------------------------------------------------------------------ the scan (plan.md)
FAMS = ("cliques", "dimer", "p4")
W_PUBLISHED = {"cliques": Fraction(1, 10), "dimer": Fraction(1, 20), "p4": Fraction(3, 10)}
K_PUBLISHED = {"cliques": (2, 3), "dimer": (2, 3), "p4": tuple(range(2, 11))}
K_TAB_FIRST = {"cliques": 3, "dimer": 3, "p4": 10}         # first row of Tab:linhas
Q_GRID = (Fraction(20, 100), Fraction(25, 100), Fraction(30, 100), Fraction(35, 100),
          Fraction(40, 100), Fraction(45, 100))
Q_STRONG = Fraction(40, 100)                                # fixed before the run
K_EXT = tuple(range(2, 11))
SCAN_MAX = 10
EXT_FAMS = ("cliques", "dimer")
G_CLOSED = {"cliques": 20.0, "dimer": 24 - 2 * math.sqrt(26), "p4": 4 - math.sqrt(5)}
CERTS = ("contour", "chain", "walk", "bridge", "CT")       # tightest first, as in the text
SRC = {"bridge": "N1", "walk": "N2", "chain": "N2p_path"}
M_CONTOUR = 2000
KAPPA_DEN = 10000                                           # kappa = j / 10^4, j = 1..9999
DPS0 = 60


def grid_w(name, q):
    """Rational bridge weight for the target q = 2w/g: exact 10q on K4 (g = 20), two decimals on
    the dimer (g = 24 - 2 sqrt 26 is irrational)."""
    if name == "cliques":
        return 10 * q
    if name == "dimer":
        return Fraction(round(float(q) * G_CLOSED["dimer"] / 2 * 100), 100)
    raise ValueError(name)


def rows():
    """(family, regime, w, q_target, k) for every row of the scan."""
    out = [(f, "published", W_PUBLISHED[f], None, k) for f in FAMS for k in K_PUBLISHED[f]]
    out += [(f, "stronger", grid_w(f, q), q, k) for f in EXT_FAMS for q in Q_GRID for k in K_EXT]
    return out


# ------------------------------------------------------------------ instances
def build_w(name, k, w):
    """families.build with every bridge weight set to w, and mu recomputed."""
    G = FAM.build(name, k)
    w = Fraction(w)
    for e in G["bridges"]:
        G["edges"][e] = w
        G["bridges"][e] = w
    G["w"] = w
    G["mu"] = max(abs(x) for x in G["edges"].values())
    return G


@contextlib.contextmanager
def _built_with(G):
    old = nm.build
    nm.build = lambda name, k: G
    try:
        yield
    finally:
        nm.build = old


def setup(name, k, w):
    """(G, P): the instance and the rigorous data of neumann.rig_setup, with checks of
    (H2), delta = w and q' = 2 max_l delta_{l,l+1} / g < 1/2 in exact rationals."""
    G = build_w(name, k, w)
    with _built_with(G):
        P = nm.rig_setup(name, k)
    assert P["delta"] == nm.fq(Fraction(w)), "delta != w"
    assert P["g"] > 2 * P["delta"], "(H2) fails"
    assert 2 * P["delta"] / P["g"] < fmpq(1, 2), "q' >= 1/2"
    return G, P


def selfcheck():
    """Verification step 1: build_w at the published weight is families.build."""
    for name in FAMS:
        for k in (2, 3, 10):
            a, b = FAM.build(name, k), build_w(name, k, W_PUBLISHED[name])
            for key in ("n", "edges", "h", "parts", "s", "t", "w", "bridges", "d", "mu"):
                assert a[key] == b[key], (name, k, key)
    for name in EXT_FAMS:
        for q in Q_GRID:
            w = grid_w(name, q)
            qa = 2 * float(w) / G_CLOSED[name]
            assert abs(qa - float(q)) < 0.002 and qa < 0.5, (name, q, w, qa)
    return True


# ------------------------------------------------------------------ certificates
def up(x):
    """Exact rational upper endpoint of an arb (or the fmpq itself)."""
    return x if isinstance(x, fmpq) else C.upper_q(x)


def closed(P, G):
    """Per-block factors {cert: {(i, x): fmpq}} of the bridge, walk-sum and chain bounds and of
    the Combes-Thomas baseline (kappa minimizing Sigma_1 on the grid), with flags."""
    out, cert = nm.bounds_123(P)
    forms, info = nm.n2prime(P)
    res, flags = {}, dict(rho_walk_certified=bool(cert["N2"]), chain_form_ok="N2p_path" in forms)
    for c, src in SRC.items():
        d = forms.get(src) if src.startswith("N2p") else out[src]
        if d is None:
            continue
        res[c] = {kk: up(v) for kk, v in d.items()}
    res["CT"], flags["CT"] = ct_best(G)
    return res, flags


def ct_best(G):
    """Combes-Thomas factors at the first kappa = j/10^4 that minimizes the upper endpoint of
    certify.sigma1 (the baseline of Section 7). A minimum over a grid fixed in advance is itself
    a bound. Returns ({(i, x): fmpq}, info)."""
    par, geo = C.parameters(G), C.geometry(G)
    his = [C.upper_q(C.sigma1(G, par, geo, fmpq(j, KAPPA_DEN))) for j in range(1, KAPPA_DEN)]
    best = min(his)
    arg = [j for j, h in enumerate(his, 1) if h == best]
    kap = fmpq(arg[0], KAPPA_DEN)
    L = ct_factors(G, par, geo, kap)
    s, t = G["s"], G["t"]
    tot = sum((L[(i, s)] * L[(i, t)] for i in range(G["k"])), fmpq(0))
    # the sum of products of per-block upper endpoints agrees with the upper endpoint of the
    # sigma1 ball to the ball radii; the verdict uses the former, itself a valid upper bound
    assert abs(float((tot - best) / best)) < 1e-20, (float(tot), float(best))
    return L, dict(kappa=str(kap), kappa_argmin_first=str(fmpq(arg[0], KAPPA_DEN)),
                   kappa_argmin_last=str(fmpq(arg[-1], KAPPA_DEN)), kappa_argmin_size=len(arg),
                   sigma1_grid_min_upper=str(best), mu=str(G["mu"]), d=G["d"],
                   eta_ct=str(C.upper_q(par["eta"])))


def ct_factors(G, par, geo, kappa):
    """The per-block factors L_i(s), L_i(t) inside certify.sigma1 (l.133-146), same formula."""
    kap = arb(kappa)
    gam = (1 + kap * par["ghat"]).log()
    c = 1 / (1 - kap)
    root = (1 - par["eta"] ** 2).sqrt()
    L = {}
    for i, V in enumerate(G["parts"]):
        A = c * par["rho"][i] * arb(len(V)).sqrt() / root
        for x, dist in ((G["s"], geo["ds"][i]), (G["t"], geo["dt"][i])):
            b = arb(1) if G["comm"][x] == i else par["eta"]
            L[(i, x)] = C.upper_q(C.amin(b, A * (-gam * dist).exp(), f"CT L_{i}"))
    return L


def contour(P, i, M=M_CONTOUR):
    """Contour certificate on block i (Proposition prop:certcontorno(b), Appendix A)."""
    r, nseg = nm.n3x_upper_reduced(P, i, M)
    return {(i, x): up(v) for x, v in r.items()}, nseg


def sigma(P, d):
    """Sigma_1 = sum_i L_i(s) L_i(t), exact rational."""
    s, t = P["s"], P["t"]
    return sum((d[(i, s)] * d[(i, t)] for i in range(P["k"])), fmpq(0))


def certifies(sig):
    """Sigma_1 < sqrt2 - 1, decided exactly: (Sigma_1 + 1)^2 < 2 for Sigma_1 >= 0."""
    return (sig + 1) ** 2 < 2


def log10q(x):
    return nm.lg10(arb(x)) if x != 0 else -math.inf


# ------------------------------------------------------------------ exact S (calibration only)
def mpq(x):
    """Fraction or fmpq -> mpf (copied from run_benchmark.mpq, l.31)."""
    if isinstance(x, fmpq):
        return mp.mpf(int(x.p)) / int(x.q)
    return mp.mpf(x.numerator) / x.denominator


def exact_once(G, P, dps):
    mp.mp.dps = dps
    H = SP.to_mp(FAM.dense(G, mpq))
    lam, Qm = SP.eig(H)
    res, orth = SP.residuals(H, lam, Qm)
    groups, min_gap, within = SP.group(lam, mp.mpf(10) ** -40)
    s, t = G["s"], G["t"]
    c = SP.coefficients(Qm, groups, s, t)
    S, Fbar = SP.S_Fbar(c)
    intervals = [(mpq(P["lo"][i] - P["delta"]), mpq(P["hi"][i] + P["delta"]))
                 for i in range(G["k"])]
    ms, cnt = SP.block_masses(Qm, lam, intervals, s)
    mt, _ = SP.block_masses(Qm, lam, intervals, t)
    tol = mp.mpf(10) ** (-(dps - 10))
    checks = dict(aglo_counts=all(cnt[i] == len(V) for i, V in enumerate(G["parts"]))
                  and sum(cnt) == G["n"],
                  residual=res <= tol, orthogonality=orth <= tol,
                  groups_separated=min_gap > mp.mpf(10) ** -30)
    return dict(S=S, Fbar=Fbar, ms=ms, mt=mt, res=res, orth=orth, min_gap=min_gap,
                within=within, n_distinct=len(groups), checks=checks)


def exact(G, P, dps0=DPS0):
    """Adaptive precision, the rule of run_benchmark.adaptive_instance (l.102-124): raise until
    S > 10^-(dps-30); if raised, confirm at dps + 20 to 1e-6 relative."""
    dps, trail = dps0, []
    while True:
        out = exact_once(G, P, dps)
        trail.append((dps, mp.nstr(out["S"], 8)))
        if out["S"] > mp.mpf(10) ** (-(dps - 30)):
            break
        est = -int(mp.floor(mp.log10(out["S"]))) if out["S"] > 0 else dps
        dps = max(dps + 60, est + 60) if out["S"] > 0 else dps * 2
    if dps > dps0:
        out2 = exact_once(G, P, dps + 20)
        rel = abs(out2["S"] - out["S"]) / out["S"]
        trail.append((dps + 20, mp.nstr(out2["S"], 8)))
        out["checks"]["precision_converged"] = bool(rel < mp.mpf(10) ** -6)
        out["S_precision_reldiff"] = mp.nstr(rel, 3)
    out["dps"], out["dps_trail"] = dps, str(trail)
    mp.mp.dps = dps
    thr = mp.sqrt(2) - 1
    out["certifies"] = bool(out["S"] < thr)
    out["margin_rel"] = mp.nstr((thr - out["S"]) / thr, 6)
    return out


# ------------------------------------------------------------------ thresholds
def thresholds(verdicts, beyond_ok=True):
    """verdicts: {k: bool} on the scanned ks; beyond_ok: the rows beyond the scan certify (the
    tabulated rows of Tab:linhas on the published families, vacuous on the extension).
    Returns (k_first, k_star, fails); k_star is None when the last scanned k fails."""
    ks = sorted(verdicts)
    fails = [k for k in ks if not verdicts[k]]
    k_first = next((k for k in ks if verdicts[k]), None)
    k_star = None
    if beyond_ok:
        for k in ks:
            if all(verdicts[kk] for kk in ks if kk >= k):
                k_star = k
                break
    return k_first, k_star, fails
