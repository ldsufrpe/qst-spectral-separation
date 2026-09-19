"""Rigorous (Arb ball arithmetic) evaluation of the Section 7 certificates.

Library: python-flint (Arb / FLINT 3).  Every arb is a midpoint-radius ball that
contains the exact value; all operations round outward.  Inputs enter as exact
rationals (fmpq) or as exact algebraic data (roots of integer polynomials enclosed
by FLINT), never as binary floats.

Formulas (manuscript/main.tex):
  eta   = delta/(g-delta)                                   (eq:eta)
  ghat  = (g-2 delta)/(2 d mu)                              (eq:eta)
  rho_i = 1 + 2(W_i+2 delta)/(pi (g-2 delta))               (eq:entrada)
  gamma = log(1+kappa ghat),  c = 1/(1-kappa)
  A_i   = c rho_i sqrt|V_i| / sqrt(1-eta^2)                 (eq:Li)
  L_i(x)= min{b_i(x), A_i exp(-gamma dist(x,V_i))},  b_i = 1 on V_i, eta off  (eq:Li)
  Sigma1= sum_i L_i(s) L_i(t)                               (eq:Ti)
  F_esp = min{1, Sigma1^2}
  b_eta = 2 eta sqrt(1-eta^2) if eta <= 1/sqrt2 else 1     (cor:comparacao)
  B_E^2 = 1 - (h_s-h_t)^2/(spread(h)+2 d mu)^2              (eq:energia)
  B^2   = min{1, Sigma1^2, b_eta^2, B_E^2}                  (eq:Bcert squared)
"""
import flint
from flint import arb, fmpq, fmpz, fmpz_mat, fmpz_poly

from families import build, dist_to_parts, bridges_form_matching, block_matrix

OVERLAPS = []          # (context) whenever a min is taken over overlapping balls


def set_prec(bits):
    flint.ctx.prec = bits


# ---------------------------------------------------------------- exact endpoints
def _arf_to_fmpq(x):
    man, exp = x.man_exp()
    man, exp = int(man), int(exp)
    return fmpq(man * 2 ** exp) if exp >= 0 else fmpq(man, 2 ** (-exp))


def upper_q(x: arb) -> fmpq:
    """Exact rational upper endpoint of the ball x (mid + rad)."""
    return _arf_to_fmpq(x.mid()) + _arf_to_fmpq(x.rad())


def lower_q(x: arb) -> fmpq:
    return _arf_to_fmpq(x.mid()) - _arf_to_fmpq(x.rad())


def amin(x: arb, y: arb, ctx=""):
    """Ball enclosing min(x, y).  Exact when the balls are ordered."""
    if upper_q(x) <= lower_q(y):
        return x
    if upper_q(y) <= lower_q(x):
        return y
    # Dominated ball (e.g. exact ties between translated blocks): if x lies
    # endpoint-wise below y, then min(x,y) is enclosed by x itself.
    if lower_q(x) <= lower_q(y) and upper_q(x) <= upper_q(y):
        return x
    if lower_q(y) <= lower_q(x) and upper_q(y) <= upper_q(x):
        return y
    OVERLAPS.append(ctx)
    lo = min(lower_q(x), lower_q(y))
    hi = min(upper_q(x), upper_q(y))
    return arb.union(arb(lo), arb(hi))


def q(x):
    """Fraction/int -> exact arb (rationals with power-of-two denominators are exact;
    others become a tight enclosing ball)."""
    return arb(fmpq(x.numerator, x.denominator)) if hasattr(x, "numerator") else arb(x)


# ------------------------------------------------------------ exact block spectra
def block_spectrum(G, i):
    """Exact char poly of H_i (integer entries in all three families), its
    factorization, and rigorous enclosures of its distinct roots."""
    M = block_matrix(G, i)
    assert all(x.denominator == 1 for row in M for x in row), "non-integer block"
    P = fmpz_mat([[int(x) for x in row] for row in M]).charpoly()
    _, factors = P.factor()
    roots = []
    for f, mult in factors:
        for r, m2 in f.complex_roots():
            assert r.imag.contains(0) and abs(r.imag).upper() < arb(2) ** -200
            roots.append((r.real, mult * m2))
    return P, factors, roots


def cluster_data(G):
    """Balls for W_i (max-min of block spectrum) and g = min_{i!=j} dist(J_i, J_j)."""
    lo, hi, facs = [], [], []
    for i in range(G["k"]):
        P, factors, roots = block_spectrum(G, i)
        facs.append(str(P) + " = " + " * ".join(f"({f})^{m}" for f, m in factors))
        vals = [r for r, _ in roots]
        mn = vals[0]
        mx = vals[0]
        for v in vals[1:]:
            mn = amin(mn, v, "block-min")
            mx = -amin(-mx, -v, "block-max")
        lo.append(mn)
        hi.append(mx)
    W = [hi[i] - lo[i] for i in range(G["k"])]
    order = sorted(range(G["k"]), key=lambda i: float(lo[i].mid()))
    # intervals disjoint and ordered => min over all pairs is min over consecutive
    gaps = [lo[order[j + 1]] - hi[order[j]] for j in range(len(order) - 1)]
    g = gaps[0]
    for x in gaps[1:]:
        g = amin(g, x, "gap-min")
    return W, g, facs


# ----------------------------------------------------------------- certificates
def parameters(G):
    assert bridges_form_matching(G)            # => ||H_out|| = w exactly (H1)
    delta = q(G["w"])
    W, g, facs = cluster_data(G)
    assert lower_q(g) > upper_q(2 * delta)     # (H2): g > 2 delta, rigorously
    d, mu = G["d"], q(G["mu"])
    eta = delta / (g - delta)
    ghat = (g - 2 * delta) / (2 * d * mu)
    rho = [1 + 2 * (Wi + 2 * delta) / (arb.pi() * (g - 2 * delta)) for Wi in W]
    return dict(delta=delta, W=W, g=g, eta=eta, ghat=ghat, rho=rho, d=d, mu=mu,
                charpolys=facs)


def geometry(G):
    ds, dxs = dist_to_parts(G, G["s"])
    dt, _ = dist_to_parts(G, G["t"])
    return dict(ds=ds, dt=dt, D=dxs[G["t"]])


def sigma1(G, par, geo, kappa: fmpq):
    kap = arb(kappa)
    gam = (1 + kap * par["ghat"]).log()
    c = 1 / (1 - kap)
    root = (1 - par["eta"] ** 2).sqrt()
    total = arb(0)
    for i, V in enumerate(G["parts"]):
        A = c * par["rho"][i] * arb(len(V)).sqrt() / root
        bs = arb(1) if G["comm"][G["s"]] == i else par["eta"]
        bt = arb(1) if G["comm"][G["t"]] == i else par["eta"]
        Ls = amin(bs, A * (-gam * geo["ds"][i]).exp(), f"L_{i}(s)")
        Lt = amin(bt, A * (-gam * geo["dt"][i]).exp(), f"L_{i}(t)")
        total += Ls * Lt
    return total


def b_eta_sq(par):
    eta = par["eta"]
    assert upper_q(eta ** 2) <= fmpq(1, 2), "eta > 1/sqrt2 branch not implemented"
    return 4 * eta ** 2 * (1 - eta ** 2)


def B_E_sq(G, par):
    h = G["h"]
    diff = q(h[G["s"]] - h[G["t"]])
    spread = q(max(h) - min(h))
    return 1 - diff ** 2 / (spread + 2 * par["d"] * par["mu"]) ** 2


def row_certificates(G, par, geo, kappa):
    S1 = sigma1(G, par, geo, kappa)
    Fesp = amin(arb(1), S1 ** 2, "F_esp")
    B2 = Fesp
    for x, tag in ((b_eta_sq(par), "b_eta"), (B_E_sq(G, par), "B_E")):
        B2 = amin(B2, x, "B2-" + tag)
    return dict(Sigma1=S1, F_esp=Fesp, B2=B2, b_eta2=b_eta_sq(par), B_E2=B_E_sq(G, par))


# ------------------------------------------------------------- 4 s.f. rounding
def _floor_log10(x: fmpq) -> int:
    assert x > 0
    e = 0
    while pow10(e) > x:
        e -= 1
    while pow10(e + 1) <= x:
        e += 1
    return e


def pow10(e):
    return fmpq(10 ** e) if e >= 0 else fmpq(1, 10 ** (-e))


def ceil_sf(x: fmpq, sf=4) -> fmpq:
    """Smallest number with `sf` significant figures that is >= x (exact)."""
    e = _floor_log10(x)
    unit = pow10(e - sf + 1)
    m = x / unit
    mi = int(m.p) // int(m.q)
    if fmpq(mi) < m:
        mi += 1
    return mi * unit


def check_printed(printed: fmpq, ball: arb, sf=4):
    """ok_upper: printed >= every value in the ball.
    ok_tight: printed is the tightest `sf`-s.f. upward rounding, i.e. the next
    smaller `sf`-s.f. number already lies below the ball's lower endpoint."""
    hi, lo = upper_q(ball), lower_q(ball)
    ok_upper = printed >= hi
    e = _floor_log10(printed)
    unit = pow10(e - sf + 1)
    # predecessor with sf figures (handles the 1.000e(e) -> 9.999e(e-1) case)
    pred = printed - unit
    if pred > 0 and _floor_log10(pred) < e:
        pred = printed - unit / 10
    ok_tight = ok_upper and pred < lo
    return ok_upper, ok_tight, ceil_sf(hi, sf)


def rounds_to(printed_str: str, ball: arb) -> bool:
    """True iff the ball lies within half a unit of the last printed digit."""
    s = printed_str.strip()
    dec = len(s.split(".")[1]) if "." in s else 0
    val = fmpq(int(s.replace(".", "").replace("-", "")), 10 ** dec)
    if s.startswith("-"):
        val = -val
    half = fmpq(1, 2 * 10 ** dec)
    return lower_q(ball) >= val - half and upper_q(ball) <= val + half
