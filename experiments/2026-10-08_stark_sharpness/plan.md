# Plan — 2026-10-08_stark_sharpness

**Date:** 2026-10-08 · **Scale:** canonical · **Modes:** numerical (exact reference, fitted rates)
and symbolic (certified interval computation)
**Parent:** opportunities.md #2 (chosen), MR1 · `derivacoes/2026-10-08_neumann_pontes/derivacao.md` §4
**Calibrated by:** `experiments/toy/2026-10-08_neumann_calibration/` (T3). That toy's plan-mini
fixed the predictions before any run, and its plan-mini and plan were approved by the author on
2026-10-08.
**Target venue:** J. Phys. A: Mathematical and Theoretical.

## Question
Is the walk-sum certificate sharp in rate on the uniform Wannier–Stark ladder? This is the one
family where a physics referee knows the exact answer: ‖P_m e_0‖ = |J_m(2J/F)|, and S ≈ sup_t|amp|
≈ (2J/F)^N/N!. The certificate is N2 and N2′, with N1, Σ₁ and N3x for comparison. The question
covers both the infinite ladder and finite chains, which are the setting of the paper's theorem.

## Hypothesis
The predictions were written before any canonical run. The toy confirmed them at F/J = 10 only.
- **Supported,** if three things hold:
  1. **On ℤ,** for F/J ∈ {8.5, 10, 16, 32, 100}, the fit log10(bound/mass) = a + b·log10 N + c·N
     over N = 10…100 gives:
     - for the per-bridge N2′, c within 2% (relative) of log10 B(q′²), with q′ = 2J/F;
     - for the ladder N2′ with global δ, c within 2% of log10(2B(q²)), with q = 4J/F;
     - b within 0.05 of ½, in both cases.

     At the S level, log10(Σ/S) − [log10(πN/2) + N·log10 B(q′²)] stays bounded (spread ≤ 0.5)
     for the per-bridge form.
  2. **The per-bridge form outside N2′'s stated hypothesis.** At F/J ∈ {5, 6} (q ≥ ½, q′ < ½),
     the per-bridge form keeps the same law. That would show the stated hypothesis q < ½ is an
     artefact of the global δ.
  3. **On finite chains** (L sites, s = 0, t = L−1, N = L−1), the N2 matrix certificate
     (Prop. N2 with δ_jl = J) is ≤ the per-bridge N2′. Its ratio to the exact finite-chain S obeys
     the same bound as on ℤ to within 1 order of magnitude, so the ends do not change the rate. The
     share of the exact rate over the last window is ≥ 0.95.
- **Refuted:**
  - log10(bound/mass) grows faster than linearly in N, so the N! is lost;
  - c misses its prediction by more than 2%;
  - the finite-chain ratio grows superlinearly.
  - **Any bound below an exact value is a correctness failure of N2′ and stops the run.**

## Variables
- **Infinite ladder (ℤ).**
  - Parameters: hₙ = nF, hopping J = 1, singleton communities, F/J ∈ {5, 6, 8.5, 10, 16, 32, 100}
    and N = 1…100. The bound is on block N with source 0.
  - The block data are exact: Wᵢ = 0, g_ij = |i−j|F, δ = ‖V‖ = 2J.
- **Finite chains.**
  - L ∈ {6, 11, 21, 41, 61} with F/J ∈ {6, 8.5, 10, 16, 32}, for 25 instances.
  - The ends are where the paper's s and t live. On a chain, ‖V‖ = 2J cos(π/(L+1)) ≤ 2J, and 2J
    is used for δ, which is safe.
- **Certificates.** Each is a per-block Lᵢ(x), with Σ = Σᵢ Lᵢ(0)Lᵢ(N):
  - **N1.**
  - **N2′ as stated:** global δ, sorted distances G₍ₘ₎ of the actual block. On ℤ these are
    ⌈m/2⌉F; on a chain they depend on the block's position.
  - **N2′, ladder form:** global δ, q^N/(2N−1)!!.
  - **N2′, per-bridge form:** coupling J on each bridge. This is the one-line variant of the proof
    that the toy README recommends for Phase 3. Here it is labelled a variant, not "N2′".
  - **N2:** the Prop. N2 matrix with δ_jl = J, on finite chains only. On singleton blocks N3
    coincides with N2, because Combes–Thomas inside a 1×1 block is the resolvent itself.
  - **N3x:** rigorous, on finite chains with L ≤ 21, M = 2000.
  - **Σ₁:** the paper's Combes–Thomas certificate, with d = 2 and μ = J for the whole chain,
    following `certify.py`, and κ optimal on a 10⁻³ grid. It is the baseline the referee called
    "extremely conservative".

  The N2′ forms are evaluated only where their walk-count series converges: q < ½ for the
  stated and ladder forms, and q′ < ½ for the per-bridge form.
- **Exact values.**
  - **On ℤ:** Bessel functions in Arb. The masses are |Jₘ(2J/F)|,
    S = Σ_{m∈ℤ}|Jₘ(x)J_{N−m}(x)|, with a rigorous tail, and sup_t|amp| = J_N(4J/F) exactly.
  - **On chains:** mpmath `eigsy` with adaptive dps, checked as follows:
    - the Lemma `aglo` count: exactly one eigenvalue in each [iF−δ, iF+δ];
    - the residual and the orthogonality;
    - S stable to 10⁻⁶ relative between dps and dps+20.
- **Dependent:**
  - log10(bound/mass) at block N and log10(Σ/S);
  - the fitted (a, b, c) on ℤ against the predictions;
  - the rate share over N-windows [10,20], [40,50] and [90,100] on ℤ, and over the L-sequence on
    chains;
  - on chains, log10(S_chain/S_ℤ), the effect of the ends.

## Method
- **Arb, 128 bits,** on ℤ. All bound formulas follow `derivacao.md` §3–§4, with ϱ̃ = 1 because
  W = 0.
- **The walk count on ℤ,** Σₘ C(N+2m, m)z^m = B(z)^N/√(1−4z), is exact for paths on ℤ. On chains
  it is an upper bound.
- **Infinite sums.** Σᵢ over ℤ is truncated at |i| ≤ N + 40, with a rigorous geometric tail bound.
- **Finite chains:**
  - the N2 matrix comes from an Arb solve of size L, with ρ < 1 certified by semipositivity;
  - N3x is computed as in `2026-10-08_neumann_14rows`, by the dense solve, since L ≤ 21;
  - the exact values come from `spectral_mp`.
- **Fits:** least squares in numpy on the Arb midpoints. The fit is descriptive only. No
  certificate depends on it.

## Methodological critique (Stage 1B)
1. **On ℤ the test checks an analytic computation.** It checks that computation honestly, but it
   is still a check of algebra. The instance with real test content is the finite chain. There
   the exact values come from diagonalization, independently of any Bessel algebra, and N2 is
   the certificate the paper can state without a variant. The README weights the finite-chain
   result accordingly.
2. **The ends of a finite chain change both S and the bound.** S_chain ≠ S_ℤ near the ends, and
   so the hypothesis is stated as a ratio to the chain's own exact S. The chain-to-ℤ ratio is
   reported to show how large the effect is.
3. **F/J = 8.5 sits near the edge q → ½.** There B(q²) ≈ 1.5 and 1/√(1−4q²) ≈ 3, so the
   global-δ forms lose most there. This is deliberate. A test only far from the edge would flatter
   the result.
4. **The fit window.** The toy gave b = 0.48 instead of ½, because the window starts at N = 5. The
   canonical fit starts at N = 10, and a tolerance of 0.05 is set in advance.

## Code reuse (read-only imports, never edited)
- The toy's `code/run.py`: `stark_bound`, `stark_S`, `stark_sigma` and `catalan_B` are copied into
  `code/stark.py`, with the SHA-256 logged and the tail bound generalized to any q < ½.
- `2026-09-18_exact_benchmarks/code/spectral_mp.py`: `to_mp`, `eig`, `residuals`, `group`,
  `coefficients`, `S_Fbar`, `block_masses`.
- `2026-10-08_neumann_14rows/code/neumann.py`: `certify_M` and the N3x routine. It is imported
  read-only after that experiment's code is frozen.
- `2026-09-18_interval_certification/code/certify.py`: the Σ₁ formulas (ρᵢ, A, γ, the η floor),
  reimplemented for singletons because `parameters()` assumes the families.

## New code
- `code/stark.py`:
  - on ℤ: the exact values and the bounds;
  - on chains: `chain_H(L, F, J)`, `chain_exact(L, F, J)`, `chain_bounds(L, F, J)` with N1, the
    N2′ forms, N2, N3x and Σ₁.
- `code/run.py`: `--part Z|chain`, `--smoke`. It appends to `experiment_log.json`.
- `code/analysis.py`: the fits, the predictions, the rate shares and the end effect, written to
  `results/{fits.csv, rate_shares.csv}`.
- `code/make_figure.py`: `figures/stark_sharpness.pdf`. Panel (a) shows log10(Σ/S) against N on ℤ
  at F/J = 10 for every certificate. Panel (b) shows the fitted c against log10 B(q′²) over F/J.

## Verification order
1. **On ℤ:** Σₙ Jₙ(x)² = 1, and S ≥ J_N(4J/F) for every N, F.
2. **On chains:** the Lemma `aglo` counts, the residual (≤ 10^−(dps−10)), the orthogonality, and
   dps convergence.
3. **Validity:** every bound is ≥ the exact value, mass and S, on ℤ and on chains.
4. **Consistency:** on chains, N2 ≤ per-bridge N2′ blockwise, and N3x ≤ N2 is expected and
   logged.
5. **Only then:** the fits, the comparison with the predictions, the shares and the figure.

## Files
- **To create:**
  - `code/{stark.py, run.py, analysis.py, make_figure.py, logutil.py}`;
  - `results/{aggregate.csv, z_series.csv, chains.csv, fits.csv, rate_shares.csv,
    viz_schema.json}`;
  - `figures/{stark_sharpness.pdf, figure_manifest.json}`;
  - `env/*`, `experiment_log.json`, `README.md`, `report_validation.md`.
- **To update:** none.
- **Read only:** as listed under Code reuse.

## Stopping criterion
- On ℤ: all 7 values of F/J, N = 1…100.
- On chains: all 25 instances. N3x on the 15 instances with L ≤ 21.

## Acceptance preview (J. Phys. A)
- One exact family, classical for physicists, where the certificate is shown sharp in rate, or
  shown not to be. Either way the answer to "extremely conservative" is concrete.
- Every bound is an Arb upper bound checked against independent exact values.

## Anticipated risks
- **Diagonalization cost at L = 61 with dps ≈ 200.** The estimate is under 2 min per instance.
- **Edge effects that make S_chain much smaller than S_ℤ.** The ratios are reported as they come
  out.
- **The per-bridge form at F/J = 5 sits where η = 2/3.** The absorption factor is then 1.34. That
  is fine, and it is reported.

## Time
- Build: about 1.5 h, after `neumann_14rows` freezes `neumann.py`.
- Run: ℤ < 1 min; chains about 30 min.
- Stage 4: about 1 h.
