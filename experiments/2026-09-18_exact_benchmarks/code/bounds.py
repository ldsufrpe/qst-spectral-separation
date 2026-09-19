"""Non-rigorous (mpmath, 50 digits) evaluation of the certificate formulas, for the
figure sweep only.  On the table rows these values are cross-checked against the
Arb balls of ../2026-09-18_interval_certification (relative agreement 1e-12).

Same formulas as certify.py of the certification experiment; see its docstring.
Block data (W_i, g, delta) are computed from the blocks, never from the spectrum
of H.
"""
import mpmath as mp

from families import bfs, block_matrix, bridges_form_matching


def block_data(G):
    """W_i, g, cluster endpoints and delta = ||H_out|| (= w for a matching)."""
    assert bridges_form_matching(G)
    ends = []
    for i in range(G["k"]):
        B = block_matrix(G, i)
        M = mp.matrix([[mp.mpf(x.numerator) / x.denominator for x in row] for row in B])
        E = mp.eigsy(M, eigvals_only=True)
        ends.append((min(E), max(E)))
    W = [hi - lo for lo, hi in ends]
    srt = sorted(ends)
    g = min(srt[j + 1][0] - srt[j][1] for j in range(len(srt) - 1))
    delta = mp.mpf(G["w"].numerator) / G["w"].denominator
    return dict(W=W, g=g, delta=delta, ends=ends)


def derived(G, bd):
    d = G["d"]
    mu = mp.mpf(G["mu"].numerator) / G["mu"].denominator
    g, delta = bd["g"], bd["delta"]
    H2 = g > 2 * delta
    eta = delta / (g - delta) if g > delta else mp.inf
    ghat = (g - 2 * delta) / (2 * d * mu)
    rho = [1 + 2 * (Wi + 2 * delta) / (mp.pi * (g - 2 * delta)) for Wi in bd["W"]]
    return dict(H2=H2, eta=eta, ghat=ghat, rho=rho, d=d, mu=mu)


def dists(G, x):
    dx = bfs(G["adj"], x)
    return [min(dx[v] for v in V) for V in G["parts"]], dx


def L_profile(G, par, kappa, x, dx_parts):
    gam = mp.log(1 + kappa * par["ghat"])
    c = 1 / (1 - kappa)
    root = mp.sqrt(1 - par["eta"] ** 2)
    out = []
    for i, V in enumerate(G["parts"]):
        A = c * par["rho"][i] * mp.sqrt(len(V)) / root
        b = mp.mpf(1) if G["comm"][x] == i else par["eta"]
        out.append(min(b, A * mp.exp(-gam * dx_parts[i])))
    return out


def certificates(G, par, kappa):
    ds, dxs = dists(G, G["s"])
    dt, _ = dists(G, G["t"])
    Ls = L_profile(G, par, kappa, G["s"], ds)
    Lt = L_profile(G, par, kappa, G["t"], dt)
    S1 = mp.fsum(a * b for a, b in zip(Ls, Lt))
    eta = par["eta"]
    b_eta = 2 * eta * mp.sqrt(1 - eta ** 2) if eta <= 1 / mp.sqrt(2) else mp.mpf(1)
    h = G["h"]
    fq = lambda x: mp.mpf(x.numerator) / x.denominator
    BE2 = 1 - (fq(h[G["s"]]) - fq(h[G["t"]])) ** 2 / (fq(max(h) - min(h)) + 2 * par["d"] * par["mu"]) ** 2
    B2 = min(1, S1 ** 2, b_eta ** 2, BE2)
    return dict(Sigma1=S1, F_esp=min(1, S1 ** 2), b_eta2=b_eta ** 2, B_E2=BE2, B2=B2,
                Ls=Ls, Lt=Lt, D=dxs[G["t"]], gamma=mp.log(1 + kappa * par["ghat"]))


def best_kappa(G, par, den=1000):
    """Grid-optimal kappa for Sigma1 (step 1/den) -- used in the sweep only."""
    # Start from kappa = 1/2 (the Section 7 convention on plateaus) and move only on a
    # strict improvement, so a flat Sigma1 returns 1/2 rather than the first grid point.
    half = mp.mpf(1) / 2
    best = (half, certificates(G, par, half)["Sigma1"])
    for j in range(1, den):
        k = mp.mpf(j) / den
        v = certificates(G, par, k)["Sigma1"]
        if v < best[1] * (1 - mp.mpf(10) ** -12):
            best = (k, v)
    return best[0]
