# Plan — 2026-09-18_binary_tree_family

**Date:** 2026-09-18 · **Scale:** canonical · **Mode:** numerical + exact combinatorics
**Parent:** opportunities.md #1 · **Delegated finding:** R-B4
**Supersedes:** `plan_mini.md` (approved 2026-09-18)

## Question
On a family whose number of communities grows *exponentially* with distance, where
Corollary `cadeiasradiais` does not apply, does Theorem C through (GR) give an
explicit prefactor `𝒞` that does not depend on k? How do `𝒞e^{−2βD}`, the sharp
`Σ1²` and the exact ITF `S²` compare?

## Construction (fixed here)
- The modules are K4 blocks (coupling 1) at the nodes of a complete binary tree of
  depth r, so k = 2^{r+1} − 1 and n = 4k. Module labels follow BFS order, `ℓ = 0…k−1`.
- **Bridges** of weight w = 0.1 join parent and child. Inside a module, vertex 0
  takes the parent bridge and vertices 1 and 2 take the child bridges, so every vertex
  carries at most one bridge. The bridges therefore form a matching, `‖H_out‖ = w`
  (H1), and d = 4 (3 internal edges plus 1 bridge), with μ = 1.
- **Site energies** are `24·ℓ` on the whole module, so `Σ_ℓ = {24ℓ−1, 24ℓ+3}`,
  g = 20 and δ = 0.1 (H2), exactly as in Example 7.1. (H3) holds with D0 = 1 and M0 = 4.
- **s** is vertex 3 of the root module, which carries no bridge. **t** is vertex 3 of
  the leftmost leaf. Then `D = 2r + 1` (1 exit hop + r bridges + (r−1) internal hops + 1 final hop; verified by BFS — the Stage 1B text said 2r+2, corrected at the smoke test).

## Hypothesis
Supported: the level step costs 2 (a bridge plus an internal edge), so
`N_x(R) ≲ C·2^{R/2}`, with `a = ln2/2 ≈ 0.347` asymptotically and C finite and
stabilizing in r. With `γ_κ = log(1 + 2.475κ)`, which is about 1.16 at κ = 0.9, there
are β with `α0 = γ_κ − β > a`, and `𝒞` is the same for every r. Meanwhile `C1^unif`,
which is linear in k, grows like 2^r.
Refuted: no admissible (κ, β, a) exists. That means the branching outruns the
certified rate, which is a finding about the reach of (GR). The fallback is a larger
Δ, which raises ĝ and γ.

## Variables
- **Combinatorial part (no diagonalization), r = 1…12:** the exact
  `C(a) := max_{x,R} N_x(R)e^{−aR}` by BFS from every vertex, for a on a grid near
  ln2/2. Also the exact `C_α = max_x Σ_i e^{−α·dist(x,V_i)}` against its (GR) bound
  `Cα/(α−a)`, then `𝒞`, `Σ1` (Arb, at the grid-optimal κ) and `C1^unif`.
- **Exact part (mpmath at 60 digits), r = 1…5 (n ≤ 252, ~7 min), plus r = 6 (n = 508,
  ~50 min) if the budget allows:** `S`, ITF `S²` and `‖P_i e_s‖` over all blocks.
- **Choice of (κ, β, a):** minimize the D at which `𝒞e^{−2βD} < 10⁻²` over a grid
  (κ ∈ {0.5,…,0.95}, β ∈ (0, γ_κ − a), a ∈ [ln2/2, 0.6]). Report the optimum and also
  a "round-number" choice for the text.

## Methodological critique of the plan-mini
1. **Theorem C is expected to turn informative only at large D.** With `c_κ = 10`,
   `A0⁴` is about 3·10⁵, so `𝒞e^{−2βD}` drops below 1 only near D ≈ 15, i.e. r ≳ 7.
   That is beyond the exact range. This is not a defect of the example: it is the
   point of (GR), because the prefactor stays fixed while k explodes. The table shows
   three layers side by side: `Σ1` (sharp, informative early), `𝒞e^{−2βD}` (uniform)
   and the exact S where computable. It must not suggest that `𝒞` is the practical
   certificate.
2. **Uniformity has to be checked, not asserted.** C(a) is computed for r up to 12
   (n ≈ 32k, BFS only), and the log shows it stabilizing.
3. **The (GR) sup over R is a sup over a finite set** (the distinct values of
   `dist(x, V_i)`), so C(a) is exact, not sampled.

## Code reuse
`families.py` and `spectral_mp.py` are copied from the sibling experiments (SHA-256
logged). New code: `tree.py` builds the family and `radial.py` does the BFS counts,
C(a), `C_α` and `𝒞`.

## Verification order
1. The matching property, d, D0 and the BFS value of D.
2. (H2), from the exact block spectra (the same K4 as in 7.1).
3. C(a) stabilizes in r.
4. `C_α,exact ≤ Cα/(α−a)`, which checks the theorem's inequality.
5. `S ≤ Σ1 ≤ C1^unif e^{−γD}` and `S² ≤ 𝒞e^{−2βD}` at every exact r. Any violation
   stops the run.
6. Double precision against mpmath where S > 1e-10.

## Files
- **To create:** `code/*`, `results/{aggregate.csv, radial_constants.csv,
  viz_schema.json}`, `experiment_log.json`, `README.md`, `env/*`.
  No figure is planned; a small table goes in §7.

## Time
Build ≈ 1.5 h. Run: BFS up to r = 12 takes minutes, and mpmath up to r = 5 takes about
10 min. r = 6 is optional, runs in the background, and takes about 50 min.
