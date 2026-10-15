"""Uniform Wannier-Stark ladder: exact values and the walk-sum certificates (see ../plan.md).

H = sum_n nF |n><n| + J sum_n (|n><n+1| + h.c.), singleton communities.
On Z (closed forms, Arb): exact masses |J_m(2J/F)|, S = sum_m |J_m(x) J_{N-m}(x)|, and
  N1, N2' as stated (sorted distances ceil(m/2)F), N2' ladder (global delta),
  N2' per-bridge (coupling J on each bridge), and the paper's Sigma1 (Combes-Thomas).
On finite chains (sites 0..L-1, s = 0, t = L-1): an instance dict compatible with
  ../../2026-10-08_neumann_14rows/code/neumann.py, so N1, N2 (Prop. N2 with
  delta_jl = J), N3 (= N2 on singletons), the N2' forms and the rigorous N3x are the
  very same code as on the 14 rows; exact values from mpmath eigsy (spectral_mp).

Ported from experiments/toy/2026-10-08_neumann_calibration/code/run.py (stark_bound,
stark_S, stark_sigma, catalan_B; SHA-256 logged), with the tail bounds generalized.
"""
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
for p in (ROOT / "experiments/2026-10-08_neumann_14rows/code",
          ROOT / "experiments/2026-09-18_exact_benchmarks/code"):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

import mpmath as mpm                                    # noqa: E402
from flint import arb, fmpq                             # noqa: E402

import neumann as nm                                    # noqa: E402
import spectral_mp as smp                               # noqa: E402

FORMS_Z = ("Sigma1", "N1", "N2p_stated", "N2p_ladder", "N2p_bridge")


# ------------------------------------------------------------------ common data
def params(F, J):
    """Exact block data: W = 0, g = F, delta = ||V|| <= 2J (2J used), rho = g/2."""
    F, J = fmpq(F), fmpq(J)
    delta = 2 * J
    assert F > 2 * delta, "(H2) fails: need F > 4J"
    return dict(F=F, J=J, delta=delta, g=F, rho=F / 2, eta=delta / (F - delta),
                q=2 * delta / F, qb=2 * J / F)


def sigma1_gamma(F, J, kappa):
    """Paper's Combes-Thomas rate and prefactor for singletons (certify.py formulas):
    ghat = (g - 2 delta)/(2 d mu), d = 2, mu = J; rho_i = 1 + 2(W + 2 delta)/(pi (g - 2 delta));
    A = c rho_i sqrt|V_i| / sqrt(1 - eta^2); L_i(x) = min{eta, A e^{-gamma |i - x|}}."""
    p = params(F, J)
    g, delta, eta = arb(p["g"]), arb(p["delta"]), arb(p["eta"])
    ghat = (g - 2 * delta) / (2 * 2 * arb(p["J"]))
    gam = (1 + arb(kappa) * ghat).log()
    rho_i = 1 + 2 * (2 * delta) / (arb.pi() * (g - 2 * delta))
    A = rho_i / (1 - arb(kappa)) / (1 - eta ** 2).sqrt()
    return A, gam


def best_kappa(F, J, N, grid=999):
    """Grid-optimal kappa (float search, step 1e-3); the bound holds for every kappa."""
    p = params(F, J)
    g, d, eta = float(p["g"]), float(p["delta"]), float(p["eta"])
    ghat = (g - 2 * d) / (4 * float(p["J"]))
    rho_i = 1 + 4 * d / (math.pi * (g - 2 * d))
    best, kb = math.inf, fmpq(1, 2)
    for a in range(1, grid + 1):
        kap = a / (grid + 1)
        A = rho_i / (1 - kap) / math.sqrt(1 - eta ** 2)
        lg_ = math.log(1 + kap * ghat)
        f = lambda n: min(eta, A * math.exp(-lg_ * n)) if n else 1.0   # noqa: E731
        tot = sum(f(abs(i)) * f(abs(N - i)) for i in range(-30, N + 31))
        if tot < best:
            best, kb = tot, fmpq(a, grid + 1)
    return kb


def z_bound(form, n, F, J, kappa=None):
    """(unfloored formula f(n), floored bound min{eta, f(n)}) for ||P_i e_s||, |i-s| = n,
    on Z; None if the form's walk-count series diverges for these F, J."""
    p = params(F, J)
    eta = arb(p["eta"])
    sq = (1 - eta ** 2).sqrt()
    q, qb = arb(p["q"]), arb(p["qb"])
    if n == 0:
        return arb(1), arb(1)
    if form == "Sigma1":
        A, gam = sigma1_gamma(F, J, kappa)
        f = A * (-gam * n).exp()
    elif form == "N1":
        f = q ** n / (1 - q) / sq
    elif form in ("N2p_stated", "N2p_ladder"):
        if not nm.lower_q(1 - 4 * q ** 2) > 0:
            return None
        pre = nm.catalan_B(q ** 2) ** n / (1 - 4 * q ** 2).sqrt() / sq
        if form == "N2p_stated":                      # G_(m) = ceil(m/2) F on Z
            prod = arb(1)
            for m in range(1, n + 1):
                prod *= q / (2 * ((m + 1) // 2) - 1)
            f = pre * prod
        else:
            f = pre * q ** n / nm.dfact(2 * n - 1)
    elif form == "N2p_bridge":
        if not nm.lower_q(1 - 4 * qb ** 2) > 0:
            return None
        f = nm.catalan_B(qb ** 2) ** n / (1 - 4 * qb ** 2).sqrt() / sq * qb ** n / nm.dfact(2 * n - 1)
    else:
        raise ValueError(form)
    return f, nm.amin(eta, f, form)


def z_ratio_bound(form, n0, F, J, kappa=None):
    """r with f(n+1) <= r f(n) for all n >= n0 (analytic, per form)."""
    p = params(F, J)
    q, qb = arb(p["q"]), arb(p["qb"])
    if form == "Sigma1":
        return (-sigma1_gamma(F, J, kappa)[1]).exp()
    if form == "N1":
        return q
    if form == "N2p_stated":
        return nm.catalan_B(q ** 2) * q / (2 * ((n0 + 2) // 2) - 1)
    if form == "N2p_ladder":
        return nm.catalan_B(q ** 2) * q / (2 * n0 + 1)
    return nm.catalan_B(qb ** 2) * qb / (2 * n0 + 1)


def z_exact_mass(n, x):
    return abs(arb(x).bessel_j(abs(n)))


def _bessel_tail(n0, h):
    """sum_{j >= n0} h^j / j!  <= 2 h^n0 / n0!   (h <= 1/2)."""
    return 2 * arb(h) ** n0 / math.factorial(n0)


def z_exact_S(N, F, J, T=40):
    """Exact S = sum_{m in Z} |J_m(x) J_{N-m}(x)|, x = 2J/F, with a rigorous tail
    (|J_n(x)| <= (x/2)^n / n!, and h = x/2 <= 1/2)."""
    x = 2 * fmpq(J) / fmpq(F)
    h = x / 2
    assert h <= fmpq(1, 2)
    tot = sum((z_exact_mass(m, x) * z_exact_mass(N - m, x) for m in range(-T, N + T + 1)), arb(0))
    tail = 2 * _bessel_tail(T + 1, h) * _bessel_tail(N + T + 1, h)
    return tot + tail * arb(0.5, 0.5)          # [tot, tot + tail]; no float conversion


def z_sigma(form, N, F, J, kappa=None, T=40):
    """Sigma = sum_{i in Z} L_i(0) L_i(N), truncated at -T <= i <= N+T, with the rigorous
    tail 2 f(T+1) f(N+T+1) / (1 - r^2)."""
    tot = arb(0)
    for i in range(-T, N + T + 1):
        a, b = z_bound(form, abs(i), F, J, kappa), z_bound(form, abs(N - i), F, J, kappa)
        if a is None or b is None:
            return None
        tot += a[1] * b[1]
    r = z_ratio_bound(form, T + 1, F, J, kappa)
    assert nm.upper_q(r) < 1
    fa, fb = z_bound(form, T + 1, F, J, kappa)[0], z_bound(form, N + T + 1, F, J, kappa)[0]
    tail = 2 * fa * fb / (1 - r ** 2)
    return tot + tail * arb(0.5, 0.5)          # [tot, tot + tail]; no float conversion


# ------------------------------------------------------------------ finite chains
def chain_P(L, F, J):
    """Instance dict in the format of neumann.rig_setup, for the chain 0..L-1."""
    p = params(F, J)
    F, J = p["F"], p["J"]
    k = n = L
    parts = [[v] for v in range(L)]
    comm = list(range(L))
    lo = hi = [v * F for v in range(L)]
    gij = [[abs(i - j) * F if i != j else None for j in range(L)] for i in range(L)]
    absV = [[J if abs(u - v) == 1 else fmpq(0) for v in range(L)] for u in range(L)]
    Vblk = [row[:] for row in absV]
    cadj = [[l for l in (j - 1, j + 1) if 0 <= l < L] for j in range(L)]
    Ncomm = [[abs(i - j) for j in range(L)] for i in range(L)]
    inner = [([[0]], 0, fmpq(0)) for _ in range(L)]
    Hblocks = [[[v * F]] for v in range(L)]
    return dict(name=f"stark_L{L}_F{F}", k=k, n=n, parts=parts, comm=comm, lo=list(lo),
                hi=list(hi), W=[fmpq(0)] * L, gij=gij, g=F, delta=p["delta"], rho=p["rho"],
                eta=p["eta"], absV=absV, Vblk=Vblk, cadj=cadj, Ncomm=Ncomm, inner=inner,
                Hblocks=Hblocks, B=list(range(L)), D=L - 1, s=0, t=L - 1)


def chain_sigma1(L, F, J, kappa):
    A, gam = sigma1_gamma(F, J, kappa)
    eta = arb(params(F, J)["eta"])
    out = {}
    for i in range(L):
        for x in (0, L - 1):
            out[(i, x)] = arb(1) if i == x else nm.amin(eta, A * (-gam * abs(i - x)).exp(), "S1")
    return out


def chain_exact(L, F, J, dps):
    """S and the block masses ||P_i e_0||, ||P_i e_{L-1}|| by mpmath eigsy at `dps`.
    Checks: one eigenvalue per I_i = [iF - delta, iF + delta] (Lemma aglo), residual,
    orthogonality."""
    mpm.mp.dps = dps
    Fm, Jm = mpm.mpf(int(F.p)) / int(F.q), mpm.mpf(int(J.p)) / int(J.q)
    H = mpm.matrix(L, L)
    for v in range(L):
        H[v, v] = v * Fm
        if v + 1 < L:
            H[v, v + 1] = H[v + 1, v] = Jm
    lam, Q = smp.eig(H)
    res, orth = smp.residuals(H, lam, Q)
    d = 2 * Jm
    owner, counts = [], [0] * L
    for j, l in enumerate(lam):
        hits = [i for i in range(L) if i * Fm - d <= l <= i * Fm + d]
        assert len(hits) == 1, "eigenvalue outside or in two intervals"
        owner.append(hits[0])
        counts[hits[0]] += 1
    assert all(c == 1 for c in counts), "Lemma aglo count fails"
    S = mpm.fsum(abs(Q[0, j] * Q[L - 1, j]) for j in range(L))
    mass = {}
    for j in range(L):
        mass[(owner[j], 0)] = abs(Q[0, j])
        mass[(owner[j], L - 1)] = abs(Q[L - 1, j])
    return dict(S=S, mass=mass, residual=res, orth=orth, lam=lam)
