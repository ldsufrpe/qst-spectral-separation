"""Round 3 (2026-09-18): prop:comparacao after the round-2 fix.

Written from scratch; imports neither the paper's code (there is none) nor
build_from_scratch.py nor the round-2 script verify_prop_comparacao.py.
The P_4-module family is rebuilt from the prose of Example ex:p4:
path on n=4k vertices, modules of 4 consecutive vertices, internal coupling 1,
bridge coupling w=0.3, site energy 5i on module i.

Claims tested (main.tex @ working tree, lines 2126-2207):
  W1  lambda_0 := max(spec H cap I_0) >= phi - delta,
      lambda_1 := min(spec H cap I_1) <= 5 - phi + delta
  A   every admissible a satisfies a <= a_star = 5/2 - phi + delta
  B   every admissible b satisfies b >= b_k = 5k - 10 + 2phi - 2delta
  E   (8.2) at (sigma, a, b) >= (8.3) for every admissible triple
  T   (8.3) > 1 for k >= 4, (8.3) at k=3 = 0.658.., > 1/2
  C   printed constants 1.444.., 0.945.., 1.074.., 0.410..
  L   (8.3)/k -> 5/(4 a_star) e^{-4 a_star/5}
  ii  |(P_0)_{ts}| <= c_kappa rho_0 e^{-gamma_kappa D}; 2.585e-9 at k=40
"""
import math
import numpy as np

PHI = (1 + math.sqrt(5)) / 2
DELTA = 0.3
A_STAR = 2.5 - PHI + DELTA


def hamiltonian(k):
    n = 4 * k
    H = np.zeros((n, n))
    for v in range(n - 1):
        H[v, v + 1] = H[v + 1, v] = 1.0 if (v + 1) % 4 else DELTA
    for v in range(n):
        H[v, v] = 5 * (v // 4)
    return H


def br_rhs(a, b, D):
    """Right-hand side of eq. (BRcomparisonbound), m = 1, |u-v| = D."""
    # q^(D/2-1/2) with q = 1 - 2a/(b+a), via log1p: plain q**n loses ~n*eps
    # relative accuracy, which swamps the O(k^-1/2) check below at k ~ 1e11.
    return 0.25 * (1 + math.sqrt(b / a)) ** 2 * math.exp(
        (D / 2 - 0.5) * math.log1p(-2 * a / (b + a)))


def b_k(k):
    return 5 * k - 10 + 2 * PHI - 2 * DELTA


def lower_83(k):
    return br_rhs(A_STAR, b_k(k), 4 * k - 1)


fails = []


def check(name, cond, detail=""):
    print(("ok   " if cond else "FAIL ") + name + ("  " + detail if detail else ""))
    if not cond:
        fails.append(name)


# ---- W1, A, B, E on the actual spectrum -----------------------------------
rng = np.random.default_rng(20260918)
worst_ratio = math.inf
max_a_seen = 0.0
min_b_minus_bk = math.inf
true_min = {}
for k in list(range(3, 41)) + [60, 100, 200, 400]:
    H = hamiltonian(k)
    lam = np.linalg.eigvalsh(H)
    lam0, lam1 = lam[3], lam[4]           # 4 lowest eigenvalues sit in I_0
    assert lam0 <= PHI + DELTA and lam1 >= 5 - PHI - DELTA  # interval sorting
    check(f"W1 k={k}", lam0 >= PHI - DELTA - 1e-12 and lam1 <= 5 - PHI + DELTA + 1e-12,
          f"lam0={lam0:.6f} lam1={lam1:.6f}")
    D = 4 * k - 1
    # admissible shifts: lam0 < sigma < lam1; extremal pair a = dist, b = max dev
    sig = np.linspace(lam0, lam1, 40003)[1:-1]
    a_ext = np.minimum(sig - lam0, lam1 - sig)
    b_ext = np.maximum(lam[-1] - sig, sig - lam[0])
    vals = np.array([br_rhs(a, b, D) for a, b in zip(a_ext, b_ext)])
    lb = lower_83(k)
    max_a_seen = max(max_a_seen, a_ext.max())
    min_b_minus_bk = min(min_b_minus_bk, (b_ext - b_k(k)).min())
    # non-extremal admissible pairs (a smaller, b larger), random
    idx = rng.integers(0, sig.size, 2000)
    a_r = a_ext[idx] * rng.uniform(0.01, 1.0, idx.size)
    b_r = b_ext[idx] * rng.uniform(1.0, 5.0, idx.size)
    vals_r = np.array([br_rhs(a, b, D) for a, b in zip(a_r, b_r)])
    m = min(vals.min(), vals_r.min())
    true_min[k] = vals.min()
    worst_ratio = min(worst_ratio, m / lb)
    check(f"E k={k}", m >= lb * (1 - 1e-12), f"min(8.2)={vals.min():.6f} (8.3)={lb:.6f}")

check("A  max admissible a <= a_star", max_a_seen <= A_STAR + 1e-12,
      f"max a={max_a_seen:.10f} a_star={A_STAR:.10f}")
check("B  admissible b >= b_k", min_b_minus_bk >= -1e-12, f"min(b-b_k)={min_b_minus_bk:.4f}")
print(f"     worst min(8.2)/(8.3) over all k and triples = {worst_ratio:.6f}")

# ---- T, C, L: closed-form claims about (8.3) -------------------------------
check("T  (8.3) at k=3 = 0.658.. > 1/2", 0.658 <= lower_83(3) < 0.659, f"{lower_83(3):.6f}")
ks = np.arange(4, 20001)
v83 = np.array([lower_83(int(k)) for k in ks])
check("T  (8.3) > 1 for 4 <= k <= 20000", v83.min() > 1, f"min={v83.min():.6f} at k={ks[v83.argmin()]}")


def expo(k):
    return (4 * k - 2) * A_STAR / (b_k(k) - A_STAR)


def g_exp(k):
    return 0.25 * (1 + math.sqrt(b_k(k) / A_STAR)) ** 2 * math.exp(-expo(k))


check("C  exponent 1.444.. at k=4", 1.444 <= expo(4) < 1.445, f"{expo(4):.6f}")
check("C  4 a_star/5 = 0.945..", 0.945 <= 4 * A_STAR / 5 < 0.946, f"{4*A_STAR/5:.6f}")
ex = np.array([expo(int(k)) for k in ks])
check("C  exponent decreasing for k>=4", np.all(np.diff(ex) < 0))
gv = np.array([g_exp(int(k)) for k in ks])
check("C  exp-minorant increasing for k>=4", np.all(np.diff(gv) > 0))
check("C  exp-minorant 1.074.. at k=4", 1.074 <= g_exp(4) < 1.075, f"{g_exp(4):.6f}")
check("C  exp-minorant <= (8.3)", np.all(gv <= v83 * (1 + 1e-12)))
lim = 5 / (4 * A_STAR) * math.exp(-4 * A_STAR / 5)
check("L  limit constant 0.410..", 0.410 <= lim < 0.411, f"{lim:.6f}")
# (1+sqrt x)^2/x = 1 + 2/sqrt x + 1/x, so the approach is O(k^{-1/2}):
# sqrt(k) * ((8.3)/k - lim) must settle at lim * 2 sqrt(a_star/5).
pred = lim * 2 * math.sqrt(A_STAR / 5)
for kk in (10**7, 10**9, 10**11):
    r = lower_83(kk) / kk
    check(f"L  (8.3)/k -> limit at O(k^-1/2), k={kk:.0e}",
          abs(math.sqrt(kk) * (r - lim) / pred - 1) < 1e-2, f"{r:.10f}")

# ---- ii: item (ii) = lema:entrada at kappa -----------------------------------
g = 5 - 2 * PHI
ghat = (g - 2 * DELTA) / (2 * 2 * 1.0)
rho0 = 1 + 2 * (2 * PHI + 2 * DELTA) / (math.pi * (g - 2 * DELTA))
gam = lambda kap: math.log(1 + kap * ghat)
check("ii gamma_1/2 = 0.1358338050", abs(gam(0.5) - 0.1358338050) < 1e-10, f"{gam(0.5):.10f}")
check("ii rho_0 = 3.0981609539", abs(rho0 - 3.0981609539) < 1e-10, f"{rho0:.10f}")
b40 = 2 * rho0 * math.exp(-159 * gam(0.5))
check("ii 2 rho_0 e^{-159 gamma} <= 2.585e-9", b40 <= 2.585e-9 and b40 > 2.58e-9, f"{b40:.6e}")
worst = 0.0
for k in (3, 5, 10, 20, 40):
    w, V = np.linalg.eigh(hamiltonian(k))
    P0 = V[:, :4] @ V[:, :4].T
    D = 4 * k - 1
    for kap in (0.1, 0.25, 0.5, 0.75, 0.9, 0.99):
        bound = rho0 * math.exp(-gam(kap) * D) / (1 - kap)
        worst = max(worst, abs(P0[D, 0]) / bound)
check("ii |(P_0)_ts| <= c_k rho_0 e^{-gamma_k D}", worst < 1, f"worst ratio={worst:.3e}")

print("\ntrue min of (8.2) over admissible triples:",
      {k: round(true_min[k], 4) for k in (3, 4, 5, 10, 40, 400)})
print("(8.3):", {k: round(lower_83(k), 4) for k in (3, 4, 5, 10, 40, 400)})
print("\nNO VIOLATION" if not fails else f"\nFAILED: {fails}")
