# Reportable narrative

On the uniform Wannier–Stark ladder (site energies $nF$, hopping $J$, one site per community), the walk-sum certificate on the graph of communities is sharp in rate once each bridge carries its own coupling. The ladder is chosen because its transfer amplitude is classical, with maximum over time about $(2J/F)^N/N!$ at distance $N$. Exact values come from multiprecision diagonalization on finite chains and from Bessel functions in ball arithmetic on the infinite ladder. There, for seven values of $F/J$ between 5 and 100, the certificate's closed form exceeds the overlap sum $S$ by a factor $(\pi N/2)\,B(q'^2)^N$, with $q'=2J/F$ and $B$ the Catalan generating function, up to an offset whose spread over $10\le N\le100$ is at most 0.004 in $\log_{10}$. The factorial decay is therefore captured. The fitted rate of its ratio to the exact block mass matches the prediction made before any run, up to a common absolute offset of $4\cdot10^{-5}$ that a $1/N$ term removes. On finite chains up to $N=60$ the certificate sits 0.9 to 2.4 orders of magnitude above $S$ and captures 0.994 to 0.996 of the exact rate, while the uniform Neumann bound and the earlier Combes–Thomas certificate sit 102 to 105 orders above at $N=60$. No certified inequality fails. Since the infinite ladder only checks an analytic computation, the finite chains carry the independent evidence. A single global coupling bound with sorted distances costs the closed form a factor $4B(q^2)$ per bridge, with $q=4J/F$, so the per-bridge statement is the sharp one.

# Operational log

**Experiment ID:** `2026-10-08_stark_sharpness`
**Date:** 2026-10-08 (start) → 2026-10-08 (end)
**Status:** complete
**Mode:** numerical (exact reference, fitted rates) + symbolic (certified interval computation)
**Opportunity:** opportunities.md #2 (chosen), MR1
**Target venue:** **J. Phys. A: Math. Theor.**, the author's decision.
**Calibrated by:** `experiments/toy/2026-10-08_neumann_calibration/` (T3). The predictions
were fixed in the toy's `plan_mini.md.archived` before any run.
**Plan:** `plan.md` (approved by the author on 2026-10-08)

## What was computed

The system is the uniform Wannier–Stark ladder H = Σₙ nF|n⟩⟨n| + J Σₙ(|n⟩⟨n+1| + h.c.), with
J = 1 and singleton communities. Its block data are exact: W = 0, g = F and δ = 2J.

- **Infinite ladder (ℤ):** F/J ∈ {5, 6, 8.5, 10, 16, 32, 100} and N = 1…100. The exact values
  are Bessel functions in Arb:
  - the masses |J_m(2J/F)|;
  - S = Σₘ|Jₘ(x)J_{N−m}(x)|, with a rigorous tail;
  - sup_t|amp| = J_N(4J/F).

  The certificates are N1, N2′ as stated, N2′ ladder (global δ), N2′ per-bridge (coupling J on
  each bridge, a one-line change in the proof of N2′) and the paper's Σ₁, with κ grid-optimal
  on a step of 10⁻³.
- **Finite chains:** sites 0…L−1, s = 0, t = L−1, N = L−1, with L ∈ {6, 11, 21, 41, 61} and
  F/J ∈ {6, 8.5, 10, 16, 32}, giving 25 instances.
  - The certificates run through the same `neumann.py` as the 14-row experiment: N1, N2
    (Prop. N2 with δ_jl = J, which on singletons equals N3), the three N2′ forms, Σ₁, and the
    rigorous N3x at M = 2000 for L ≤ 21.
  - The exact values come from mpmath `eigsy` with adaptive dps (70–214), checked against
    dps + 20, the Lemma `aglo` count, the residual and orthogonality.

## Numerical results

**Infinite ladder.** The fit is log10(bound/exact mass) = a + b·log10 N + c·N over
N = 10…100, against the predictions fixed in advance. The S-level offset and spread are those
of log10(Σ/S) − [log10(πN/2) + N·log10 B(q′²)] for the per-bridge form.

| F/J | per-bridge b | per-bridge c (pred.) | ladder c (pred.) | as-stated c (pred.†) | S-level offset / spread |
|---|---|---|---|---|---|
| 5 | 0.490 | 0.09696 (0.09691) | — | — | +0.701 / 0.0040 |
| 6 | 0.491 | 0.05920 (0.05915) | — | — | +0.381 / 0.0030 |
| 8.5 | 0.491 | 0.02637 (0.02633) | 0.47568 (0.47563) | 0.77669 (0.77666) | +0.153 / 0.0019 |
| 10 | 0.492 | 0.01856 (0.01852) | 0.39798 (0.39794) | 0.69900 (0.69897) | +0.104 / 0.0014 |
| 16 | 0.492 | 0.00699 (0.00695) | 0.33119 (0.33114) | 0.63220 (0.63217) | +0.037 / 0.0005 |
| 32 | 0.492 | 0.00175 (0.00171) | 0.30802 (0.30798) | 0.60904 (0.60901) | +0.009 / 0.0001 |
| 100 | 0.492 | 0.00022 (0.00017) | 0.30177 (0.30173) | 0.60278 (0.60276) | +0.001 / 0.0000 |

† The as-stated prediction log10(4B(q²)) was derived after the toy, from its F/J = 10 value.
The per-bridge and ladder predictions predate any run. For the ladder and as-stated forms,
b is 0.491–0.497. At F/J = 5 and 6 those two forms fall outside their hypothesis q = 4J/F < ½,
so they are not evaluated.

**Distance to S at N = 100**, log10(Σ/S):

| F/J | Σ₁ | N1 | N2′ as stated | N2′ ladder | N2′ per-bridge |
|---|---|---|---|---|---|
| 5 | 194.1 | 191.7 | — | — | 12.6 |
| 10 | 194.2 | 190.5 | 72.3 | 42.5 | 4.2 |
| 100 | 194.3 | 190.1 | 62.2 | 32.4 | 2.2 |

**Finite chains,** log10(Σ/S_chain). The prediction column is the infinite-ladder per-bridge
law log10(πN/2) + N·log10 B(q′²).

| L (N) | F/J | end effect log10(S_chain/S_ℤ) | N2 | prediction | N3x | N1 | Σ₁ |
|---|---|---|---|---|---|---|---|
| 6 (5) | 6 / 10 / 32 | −0.005 / −0.002 / −0.000 | 1.29 / 1.01 / 0.91 | 1.19 / 0.99 / 0.90 | 0.43 / 0.24 / 0.17 | 4.8 / 4.6 / 4.4 | 4.8 / 5.3 / 5.7 |
| 21 (20) | 6 / 10 / 32 | −0.000 | 1.90 / 1.62 / 1.51 | 2.68 / 1.87 / 1.53 | 0.74 / 0.54 / 0.47 | 26.7 / 26.2 / 25.8 | 28.0 / 28.3 / 28.5 |
| 61 (60) | 6 / 10 / 32 | −0.000 | 2.38 / 2.09 / 1.98 | 5.52 / 3.09 / 2.08 | — | 102.8 / 102.2 / 101.9 | 105.3 / 105.4 / 105.5 |

All 25 instances are in `results/chain_summary.csv` and `results/chains.csv`.

**Share of the exact rate.**

| certificate | ℤ, last window N = 90→100, 7 values of F/J | chains, L = 21→61, 5 values of F/J |
|---|---|---|
| Σ₁ | 0.03–0.38 | 0.06–0.31 |
| N1 | 0.04–0.38 | 0.08–0.32 |
| N2′ per-bridge (ℤ) / N2′ path (chains) | 0.957–0.999 | 0.966–0.995 |
| N2 (chains) | — | 0.994–0.996 |

**Checks:**
- **Violations: 0.** On ℤ, every form is ≥ the exact mass at every distance 0…141 and Σ ≥ S
  for every N. On chains, every certificate is ≥ the exact mass in every block, and Σ ≥ S, on
  all 25 instances.
- **Exact values on ℤ:** Σₙ Jₙ(x)² = 1 holds to the working precision, and S ≥ J_N(4J/F) holds
  rigorously for every N and F/J.
- **Exact values on chains:** the Lemma `aglo` count is exactly one eigenvalue per
  [iF−δ, iF+δ]. The residual is ≤ 10^−(dps−10) and the S relative change at dps + 20 is
  < 10⁻⁶ on all 25.
- **Consistency on chains:** N2 ≤ N2′ path blockwise holds on all 25, and ρ < 1 is certified
  for every walk matrix.
- **Cost:**
  - ℤ part: about 4 min serial;
  - chains: about 8 min on 6 workers. L = 61 diagonalization takes about 63 s, and N3x at
    L = 21 takes about 190 s.

## Hypothesis check (`plan.md`)

1. **Infinite ladder: supported, with one criterion failing as written.**
   - b lies within 0.05 of ½ for every form and F/J, at 0.490–0.492.
   - The ladder form's c lies within a relative 0.02 of log10(2B(q²)) at every F/J ≥ 8.5. The deviations
     are at most +5·10⁻⁵ absolute, about 10⁻⁴ relative.
   - The per-bridge form's c lies within a relative 0.02 of log10 B(q′²) at five of the seven
     F/J. It **fails the relative-0.02 criterion at F/J = 32 (+0.025) and F/J = 100 (+0.24).**
   - The cause is a common absolute offset of +4·10⁻⁵ in c, the same in every form and every
     F/J, which is an O(1/N) correction the plan's fit form omits. Where the predicted c is
     only 1.7·10⁻³ and 1.7·10⁻⁴, that offset is large in relative terms.
   - Adding a 1/N term to the fit brings c to within 10⁻⁶ absolute of its prediction for
     the per-bridge and ladder forms at every F/J (`fits.csv: c_with_1overN`). The as-stated
     form then sits −8.8·10⁻⁵ off its post-toy prediction. Its polynomial prefactor differs, and
     that prediction was not fixed in advance.
   - The S-level law holds with spread ≤ 0.004 over N = 10…100 at every F/J.
2. **The per-bridge form at F/J ∈ {5, 6}: supported.** The same law holds where the stated N2′
   does not apply (q ≥ ½). The hypothesis q < ½ of N2′ is therefore an artefact of the global δ.
3. **Finite chains: supported.**
   - N2 sits below its prediction + 1 on all 25 instances, and in fact below the prediction
     itself from L = 11 on.
   - The ends do not move S (|log10(S_chain/S_ℤ)| ≤ 0.005).
   - N2 captures 0.994–0.996 of the exact rate between L = 21 and L = 61.

**MR1 answered.** On the canonical physics family, the walk-sum certificate is sharp in rate.
- At N = 60 it sits 2.0–2.4 orders above the exact S of the finite chain, while N1 and the
  paper's Σ₁ sit 102–105 orders above.
- On ℤ the per-bridge closed form reproduces (2J/F)^N/N! up to (πN/2)·B(q′²)^N, as the napkin
  estimate of 2026-10-08 said.
- N2′ exactly as written in `derivacao.md` (global δ, sorted distances) loses a factor
  4B(q²) per bridge. That is 72 orders at N = 100 and F/J = 10. Stating N2′ with per-bridge
  couplings removes this loss (a Phase 3 decision, not acted on here).

**Caveat.** On ℤ this checks an analytic computation on a closed-form family. The independent
content is the finite-chain result, whose exact values come from diagonalization.

## Figures

- `figures/2026-10-08_stark_sharpness_rate.pdf`, **for the manuscript** (with PNG and
  `figure_manifest.json`).
  - log10(Σ/S) against N at F/J = 10: every certificate on ℤ, the asymptotic law (with its
    constant 1/((1−4q′²)(1−η²)) since S9), and N2 and N3x on the finite chains.
  - **Single panel since Phase 3 S9 (2026-10-13), by the author's choice.** The former panel (b),
    fitted against predicted c over all F/J, was dropped. The manuscript states the law at every
    F/J through its maximum deviation (`manuscript/make_experimental_macros.py`, StarkLawDevMax).
  - **Legend labels set in S9** to the manuscript's names (`LABELS` in `code/make_figure.py`):
    N2′ per-bridge → chain bound (ii); N2′ ladder, global δ → chain bound (ii), δ_jl = 2J; N2′ as
    stated → chain bound (i), δ_jl = 2J; N2 → walk-sum bound; N3x → contour certificate;
    N1 → bridge bound; Σ₁ → Combes–Thomas bound. It is Figure `F:stark` of the manuscript.
  - The width is 6.0 in. The IOP class is not installed, so the J. Phys. A text width was not
    measured.

## Pointers

- Aggregate data: `results/aggregate.csv` (long format: part, F/J, N, certificate, log10 Σ,
  log10 S, distance), `results/z_series.csv`, `results/chains.csv`, `results/chain_blocks.csv`,
  `results/fits.csv`, `results/rate_shares.csv` and `results/chain_summary.csv`.
- Superseded: `results/z_series_pass1.csv` and `results/run_Z_pass1.log` (see Notes).
- Full log: `experiment_log.json`, 43 runs.
  - r0001–r0004: smoke.
  - r0005–r0011: ℤ pass 1, superseded. r0010 and r0011 have status failure.
  - r0012–r0018: ℤ final.
  - r0019–r0043: chains.
- Code: `code/stark.py`, `code/run.py`, `code/analysis.py` and `code/make_figure.py`. It reuses
  `../2026-10-08_neumann_14rows/code/neumann.py` and `../2026-09-18_exact_benchmarks/code/spectral_mp.py`
  read-only, with SHA-256 logged.
- Reproduce:
  ```bash
  cd experiments/2026-10-08_stark_sharpness/code
  ../../../.venv/bin/python -I run.py --part Z       # ~4 min
  ../../../.venv/bin/python -I run.py --part chain   # ~8 min, 6 workers
  ../../../.venv/bin/python -I analysis.py && ../../../.venv/bin/python -I make_figure.py
  ```

## Notes

- **Smoke r0002 failed on a real bug, and the validity check caught it.** In the reduced N3x
  system of `neumann.py`, the |V|U row of a vertex carrying two bridges was overwritten. That
  produced zero bounds, below the exact masses. The fix accumulates over bridge partners, and
  the crosscheck of the 14-row experiment now includes a non-matching instance.
- **ℤ pass 1 (r0005–r0011) had a floor of 10⁻³⁰⁰,** added to the tail radius to avoid a float
  underflow. Values of S below about 10⁻³⁰⁰ (F/J = 100, N ≥ 92) were pinned to that floor. The
  S ≥ sup|amp| check failed on F/J = 32 and 100 and exposed it. The tail is now added as the
  ball tail·[0, 1], with no float conversion, and r0012–r0018 are the final runs. The toy carried
  the same floor, but all its values lie above 10⁻⁸⁸, so T3 is unaffected.
- **Σ₁ on singletons** uses the paper's formulas from `certify.py` with d = 2 and μ = J for the
  whole chain. Its κ is chosen in float on a grid of step 10⁻³, and the bound is then
  evaluated in Arb, so it is valid for every κ.
