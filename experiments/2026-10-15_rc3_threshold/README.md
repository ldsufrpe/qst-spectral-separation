# Reportable narrative

The refined certificates were expected to guarantee an end-to-end transfer fidelity below the classical value 2/3 on shorter modular chains than the Combes–Thomas baseline. Every bound was evaluated in ball arithmetic from block data and its verdict Σ₁ < √2 − 1 decided in rational arithmetic, alongside the exact overlap sum S, a multiprecision calibration, on the three published families and on K4 and dimer chains with stronger bridges, q = 0.20–0.45 and k = 2–10. On the published P4 chains, S and the contour, walk-sum and chain bounds fall below √2 − 1 from k = 2, the bridge bound from 4 and the baseline only from 9, failing at k = 3–8. At q = 0.40 the contour, like S, certifies from k = 2, the other refined bounds need k = 3 to 5, and the baseline needs 7 on the dimer and more than 10 on K4. Every certificate bounds S from above on all 121 chains. The refined certificates thus change the baseline's verdict, which the published K4 and dimer chains cannot show because η ≈ 0.004–0.005 caps each block factor and every bound certifies at k = 2. On their stronger bridges the change first appears at q = 0.30, where the baseline's first failure exceeds √2 − 1 by about 10⁻³. The evidence rests instead on P4 and on q = 0.40, where μ = w also penalizes the baseline.

# Operational log

**Experiment ID:** `2026-10-15_rc3_threshold`
**Date:** 2026-10-15 (start) → 2026-10-15 (end)
**Status:** complete
**Mode:** numerical (exact S in multiprecision) + symbolic (certificates in Arb, verdicts in exact
rationals)
**Opportunity:** opportunities.md #2; review card **R-C3** of Phase 6, delegated on 2026-10-15
**Target venue:** **J. Phys. A: Math. Theor.**, the author's decision.
**Plan:** `plan.md`, with Stage 1A in `plan_mini.md.archived`. The author approved both on
2026-10-15, choosing the stable threshold k\*, the extension on K4 and the dimer only, six
values of q and no figure.

## What was computed

For each chain, s and t sit at the two ends, as in §7. Each certificate is Σ₁ of
`teo:principal`, with per-block factors that are upper endpoints of Arb balls, held as exact
rationals. Its verdict is (Σ₁ + 1)² < 2, i.e. Σ₁ < √2 − 1, so F_Bloch < 2/3 at every time by
`cor:fent`. The test is exact.
- **Certificates.** The code of the published campaigns is imported read-only, and no file in
  those folders was edited:
  - bridge (`teo:pontes`, N1);
  - walk-sum (`prop:passeios`, N2);
  - chain (`eq:cadeiamassa` with the coupling of each bridge, N2p_path);
  - contour (`prop:certcontorno`(b), M = 2000 arcs, 128 bits);
  - Combes–Thomas (the baseline of §7, κ minimizing Σ₁ on j/10⁴).
- **Exact S** comes from mpmath `eigsy` at 60 to 120 digits, by the adaptive rule of
  `2026-09-18_exact_benchmarks`. It is calibration only, with numerical verdicts. The largest S
  over the scan is 7.6·10⁻³, below 0.02 of the threshold, so no verdict on S is close.
- **Rows.** There are 121 of them:
  - **published** (`Tab:familias`): K4 k = 2, 3; dimer k = 2, 3; P4 k = 2…10. The upper end is
    the first row of `Tab:linhas`;
  - **stronger bridges** on K4 and the dimer, k = 2…10. Only w changes, on the grid
    q = 2w/g ∈ {0.20, 0.25, 0.30, 0.35, 0.40, 0.45}:
    - K4: w = 2 … 4.5;
    - dimer: w = 1.38 … 3.11, two decimals.

    Every row satisfies g > 2δ and q′ = q < 1/2. Here w > 1, so μ = w.
- **Thresholds.** For each family, w and certificate:
  - k_first, the first k that certifies;
  - **k\***, the smallest k from which it certifies at every scanned k. On the published
    families, the rows of `Tab:linhas` beyond the scan count too (K4 4–10, dimer 5–12, P4 20 and
    40), and they all certify;
  - the k that fail.

  The **onset** is the smallest q of the grid with a k ≤ 10 at which the baseline fails and the
  contour certifies. It was fixed as a rule before the run, as was the comparison point
  q = 0.40.

## Numerical results

**P4, published parameters (η = 0.205, μ = 1).** Each cell is Σ₁ (the exact S in the first
column), with ✓ when it lies below √2 − 1 ≈ 0.414.

| k | exact S | contour | chain | walk-sum | bridge | CT |
|---|---|---|---|---|---|---|
| 2 | 0.002 ✓ | 0.005 ✓ | 0.410 ✓ | 0.410 ✓ | 0.410 ✓ | 0.410 ✓ |
| 3 | 5·10⁻⁷ ✓ | 8·10⁻⁶ ✓ | 0.182 ✓ | 0.131 ✓ | 0.452 ✗ | 0.452 ✗ |
| 4 | 3·10⁻¹¹ ✓ | 1·10⁻⁹ ✓ | 0.033 ✓ | 0.023 ✓ | 0.348 ✓ | 0.494 ✗ |
| 5 | 6·10⁻¹⁶ ✓ | 6·10⁻¹⁴ ✓ | 0.006 ✓ | 0.003 ✓ | 0.186 ✓ | 0.536 ✗ |
| 6 | 5·10⁻²¹ ✓ | 1·10⁻¹⁸ ✓ | 3·10⁻⁴ ✓ | 2·10⁻⁴ ✓ | 0.103 ✓ | 0.578 ✗ |
| 7 | 2·10⁻²⁶ ✓ | 9·10⁻²⁴ ✓ | 1·10⁻⁵ ✓ | 5·10⁻⁶ ✓ | 0.053 ✓ | 0.620 ✗ |
| 8 | 6·10⁻³² ✓ | 4·10⁻²⁹ ✓ | 3·10⁻⁷ ✓ | 1·10⁻⁷ ✓ | 0.024 ✓ | 0.425 ✗ |
| 9 | 1·10⁻³⁷ ✓ | 1·10⁻³⁴ ✓ | 7·10⁻⁹ ✓ | 2·10⁻⁹ ✓ | 0.010 ✓ | 0.317 ✓ |
| 10 | 1·10⁻⁴³ ✓ | 1·10⁻⁴⁰ ✓ | 1·10⁻¹⁰ ✓ | 4·10⁻¹¹ ✓ | 0.004 ✓ | 0.247 ✓ |

**Stable threshold k\*** (`results/thresholds.csv`):

| family, regime | q | exact S | contour | chain | walk-sum | bridge | CT |
|---|---|---|---|---|---|---|---|
| K4, published | 0.010 | 2 | 2 | 2 | 2 | 2 | 2 |
| dimer, published | 0.007 | 2 | 2 | 2 | 2 | 2 | 2 |
| P4, published | 0.340 | 2 | 2 | 2 | 2 | 4 | 9 |
| K4, q = 0.20 and 0.25 | | 2 | 2 | 2 | 2 | 2 | 2 |
| K4, q = 0.30 (onset) | 0.300 | 2 | 2 | 2 | 2 | 2 | 8 |
| K4, q = 0.35 | 0.350 | 2 | 2 | 3 | 3 | 4 | 10 |
| K4, q = 0.40 | 0.400 | 2 | 2 | 3 | 3 | 4 | > 10 |
| K4, q = 0.45 | 0.450 | 2 | 2 | 4 | 3 | 5 | > 10 |
| dimer, q = 0.20 and 0.25 | | 2 | 2 | 2 | 2 | 2 | 2 |
| dimer, q = 0.30 (onset) | 0.300 | 2 | 2 | 2 | 2 | 2 | 5 |
| dimer, q = 0.35 | 0.351 | 2 | 2 | 3 | 3 | 4 | 6 |
| dimer, q = 0.40 | 0.400 | 2 | 2 | 3 | 3 | 5 | 7 |
| dimer, q = 0.45 | 0.451 | 2 | 2 | 4 | 3 | 6 | 10 |

The k at which the baseline fails while the contour certifies:
- P4 published: 3–8;
- K4: 4–7 at the onset and 2–10 at q = 0.40;
- dimer: 4 at the onset and 2–6 at q = 0.40.

**Checks:**
- **Validity:** 0 violations in 121 rows × 5 certificates. Σ₁ ≥ S always, and every per-block
  factor is ≥ its exact mass (to 10⁻⁶ relative).
  - The `aglo` counts, the residual and the orthogonality pass on every row.
  - ρ < 1 is certified for every walk matrix and every contour arc.
- **Reproduction:** 18/18 on the rows that `Tab:linhas` already has (K4 k = 3, dimer k = 3,
  P4 k = 10), in `results/reproduction.csv`:
  - bridge, walk-sum, chain and contour are identical to `2026-10-08_neumann_14rows`;
  - CT is within 3·10⁻⁷ in log10 of `2026-09-18_interval_certification`;
  - S is within 2·10⁻⁷ relative of `2026-09-18_exact_benchmarks`.
- **Probe:** the verdicts agree with the pre-plan probe on 294 of 294 shared (row, certificate)
  pairs (`results/probe_agreement.csv`).

## Hypothesis check (`plan.md`)

- **P4 published:** supported exactly. k\* is 2 for the contour, chain, walk-sum and S; 4 for the
  bridge bound (fails at k = 3); 9 for the baseline (fails at k = 3…8). The refined certificates
  change the baseline's verdict on six chain lengths of a published family, at μ = 1.
- **K4 and dimer:** supported.
  - The onset is q = 0.30 on both.
  - From q = 0.35 on, 2η > √2 − 1, and at k = 2 only the contour certifies, as S does.
  - At q = 0.40 the baseline needs k = 7 on the dimer and more than 10 on K4. The bridge bound
    needs 4–5 and the walk-sum and chain bounds need 3.
- **Refutation criteria:** none met. The delimitation route of R-C3 is not needed.

## Figures

None, by the decision of Stage 1A. The two tables above carry the answer.
`results/viz_schema.json` gives the dashboard three views of `aggregate.csv` (P4 published, K4
and dimer at q = 0.40). They are not manuscript figures.

## Pointers

- **Code:**
  - `code/rc3.py`: instances, certificates, the exact S and the thresholds;
  - `code/run.py`: the pool and the log;
  - `code/aggregate.py`: CSVs, checks, thresholds and the summary;
  - `code/ct_refine.py`;
  - `code/logutil.py`.
- **Results:**
  - `results/aggregate.csv`, one row per (family, w, k, certificate);
  - `results/thresholds.csv`;
  - `results/reproduction.csv`;
  - `results/probe_agreement.csv`;
  - `results/ct_refine.csv`;
  - `results/raw_final.json`, the exact rational per-block bounds;
  - `results/run_final.log` and `results/run_repro.log`.
- **Log:** `experiment_log.json`, 123 runs, none failed:
  - r0001 smoke;
  - r0002–r0122, one per row;
  - r0123 ct_refine.

  Its `provenance_note` cites R-C3.
- **Environment:** `env/`. Python 3, python-flint 0.9.0 (Arb) and mpmath 1.3.0, on 6 workers.
- **Macros:** `manuscript/make_experimental_macros.py` reads `thresholds.csv` and `aggregate.csv`
  and writes the `\Threshold*` macros of `manuscript/experimental_macros.tex`. Generating them
  changed no other macro (`--check-unchanged`). `main.tex` uses none yet, and its PDF text is
  identical before and after.
- **Pre-plan probe:** `exploration/preplan_probe/`. It used M = 300, κ step 1/50 and float64 S.
  It shaped the design and is not evidence.

## Notes

- **The η floor.** Every certificate takes min{η, ·} per block, so at k = 2 it gives Σ₁ ≤ 2η.
  - On the published K4 and dimer, η ≤ 0.005, so every certificate certifies at every k. The
    largest Σ₁ there is 0.010.
  - On P4, 2η = 0.410, only 0.004 below 0.414. Everything certifies at k = 2 and the bridge and
    baseline fail at k = 3, which is why the threshold is k\* (stable) rather than the first k.
- **At k = 2, the bridge, walk-sum and chain bounds and the baseline all equal 2η.** With one
  bridge, their factor exceeds η. Only the contour, which reads the entries of the block resolvents, goes
  below η. Where 2η > √2 − 1, k = 2 therefore separates the contour from the rest by
  construction.
- **The onset sits on the η floor of the baseline.** At q = 0.30 the baseline first fails at
  k = 4 on both families, with Σ₁ = 0.415, about 10⁻³ above the threshold. On K4 it keeps
  failing up to k = 7.
  - On those first failures (and at k = 2 for q = 0.35), its factor is capped at η for every κ,
    so Σ₁ = 2η + (k − 2)η², independent of κ.
  - `code/ct_refine.py` re-evaluates the five CT failures below 1.05 times the threshold on
    κ = j/10⁶. All still fail, and four are η-capped for every κ. This check was added after
    the run.
  - For the paper, P4 at the published parameters and q = 0.40 are the cases that do not hinge
    on a near-tie. There the failing baseline lies at 1.03 to 2.26 times the threshold.
- **μ = w on the extension.** Lemma `lema:ct` needs μ ≥ every coupling, bridges included, so a
  strong bridge also lowers the baseline's rate per edge. This is the lemma as stated, and it
  confounds the extension. P4 at the published parameters (μ = 1) shows the verdict change
  without it.
- **S versus the supremum of the amplitude.** S ≥ sup|f|, so the exact column bounds the true
  crossing from above. It is 2 everywhere, the smallest k with bridges, so the estimate is exact.
- **k\* is relative to the scan.** Nothing is claimed beyond k = 10 on the extension.
- **The "Pendência da Phase 4" was corrected in this session,** in documentation only:
  - `2026-09-18_exact_benchmarks/plan.md`, the F_samp item;
  - `2026-10-08_stark_sharpness/README.md`, the dps comparison.
