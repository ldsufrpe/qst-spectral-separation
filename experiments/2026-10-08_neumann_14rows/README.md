# Reportable narrative

Upper bounds on the transfer quantity S(s,t) = Σ_λ|⟨s|Π_λ|t⟩| built from block-local data alone can capture nearly all of its exponential decay in the number k of communities and lie far closer to S than the earlier Combes–Thomas certificate. On fourteen chains of K4 modules, dimers or P4 modules, every bound is computed in ball arithmetic from block spectra, block resolvents and bridge weights, and exact values of S from multiprecision diagonalization serve only for calibration. The contour integral with exact block resolvents, evaluated as an upper Riemann sum, lies 0.6 to 5.5 orders of magnitude above S and at most 2.0 above the best that any bound built from block masses can reach. The earlier certificate lies 5.6 to 259.5 orders above S, and the best closed form, a hybrid bound with a Neumann series across blocks and Combes–Thomas estimates inside them, 3.3 to 80.7. The contour integral captures 0.98 to 0.99 of the exact decay rate in k, the hybrid bound 0.72 to 0.83 and the earlier certificate 0.05 to 0.17. No certified inequality fails on any instance. The hybrid bound therefore recovers most of the exact decay rate but stays far above S. Only the contour integral is close in both respects, and it is a computer-assisted certificate evaluated instance by instance, not an a priori formula.

# Operational log

**Experiment ID:** `2026-10-08_neumann_14rows`
**Date:** 2026-10-08 (start) → 2026-10-08 (end)
**Status:** complete
**Mode:** numerical (benchmark against exact values) + symbolic (certified interval computation)
**Opportunity:** opportunities.md #2 (chosen), MR2 and MR3
**Target venue:** **J. Phys. A: Math. Theor.**, the author's decision.
**Calibrated by:** `experiments/toy/2026-10-08_neumann_calibration/` (T1, T2)
**Plan:** `plan.md` (approved by the author on 2026-10-08)

## What was computed

Every certificate below is a rigorous Arb upper bound on S(s,t) = Σ_λ|⟨s|Π_λ|t⟩|. The inputs
are exact:
- the block hulls Jᵢ come from the exact characteristic polynomials;
- g, g_ij, Wᵢ and η are rationals on the safe side;
- δ = w exactly, because the bridges form a matching (asserted).

It was evaluated on the 14 rows of Table 3, with s and t as in §7.

- **Σ₁:** the paper's Combes–Thomas certificate, read from
  `2026-09-18_interval_certification`.
- **N1, N2, N2′ and N3:** the closed forms of `derivacoes/2026-10-08_neumann_pontes/derivacao.md`.
  - N2′ is reported in three forms: as stated (sorted distances), path (the proof's own
    first-arrival product) and ladder. The ladder condition was certified on all 14 rows.
  - N3 is minimized over κ ∈ {1/2, 7/10, 4/5, 9/10, 19/20}, a grid fixed in the plan.
- **N3x:** a rigorous upper Riemann sum over each contour Γᵢ, at M = 1000 and 2000 segments.
  - Uses only block data: block resolvents along Γᵢ, δ and the bridges.
  - **It is a per-instance computer-assisted certificate, not a closed-form bound.**
- **best:** the per-block minimum over all certificates. Since each one is a valid bound, so is
  their minimum.

The exact S and block masses come from `2026-09-18_exact_benchmarks`. They are used for
calibration only, and no certificate reads them.

## Numerical results

**Distance to the exact value, log10(Σ/S)** (`results/aggregate.csv`). The floor is
Σᵢ‖Pᵢeₛ‖‖Pᵢeₜ‖ with exact masses, the best any block-mass certificate can do.

| row | D | log10 S | Σ₁ | N1 | N2 | N2′ stated | N2′ path | N2′ ladder | N3 | **N3x** | floor |
|---|---|---|---|---|---|---|---|---|---|---|---|
| K4 k=3 | 5 | −7.6 | 5.6 | 4.0 | 3.6 | 3.6 | 3.6 | 3.6 | 3.3 | **0.6** | 0.3 |
| K4 k=4 | 7 | −11.8 | 9.8 | 6.4 | 5.5 | 5.9 | 5.5 | 5.5 | 4.9 | **0.7** | 0.3 |
| K4 k=5 | 9 | −16.2 | 13.8 | 8.9 | 7.5 | 8.5 | 7.5 | 7.6 | 6.5 | **0.8** | 0.3 |
| K4 k=6 | 11 | −20.8 | 17.5 | 11.6 | 9.5 | 10.7 | 9.5 | 9.7 | 7.9 | **0.9** | 0.3 |
| K4 k=8 | 15 | −30.5 | 25.8 | 17.4 | 13.6 | 15.5 | 13.6 | 13.9 | 10.3 | **1.1** | 0.4 |
| K4 k=10 | 19 | −40.6 | 34.5 | 23.6 | 17.8 | 20.3 | 17.8 | 18.2 | 12.7 | **1.3** | 0.4 |
| dimer k=3 | 5 | −8.5 | 6.4 | 4.8 | 4.2 | 4.2 | 4.2 | 4.4 | 3.7 | **1.2** | 1.0 |
| dimer k=5 | 9 | −17.7 | 14.5 | 10.0 | 8.4 | 9.6 | 8.4 | 8.8 | 5.9 | **1.5** | 1.2 |
| dimer k=7 | 13 | −27.7 | 22.9 | 15.9 | 12.5 | 14.4 | 12.5 | 13.3 | 7.7 | **1.8** | 1.4 |
| dimer k=9 | 17 | −38.2 | 31.5 | 22.2 | 16.7 | 19.2 | 16.7 | 17.9 | 9.2 | **1.9** | 1.5 |
| dimer k=12 | 23 | −54.6 | 44.6 | 32.4 | 23.2 | 26.7 | 23.2 | 25.0 | 11.4 | **2.1** | 1.6 |
| P4 k=10 | 39 | −43.0 | 42.4 | 40.6 | 32.6 | 36.4 | 33.2 | 35.9 | 22.0 | **3.2** | 1.9 |
| P4 k=20 | 79 | −110.7 | 107.6 | 104.1 | 80.1 | 87.9 | 81.3 | 88.4 | 41.9 | **4.3** | 2.7 |
| P4 k=40 | 159 | −270.4 | 259.5 | 254.7 | 192.6 | 207.9 | 195.1 | 211.0 | 80.7 | **5.5** | 3.5 |

**Share of the exact decay rate in k** (MR2, `results/rate_shares.csv`). The share is
slope(log10 Σ)/slope(log10 S), in per cent, computed three ways: endpoint secant /
least squares / last pair.

| certificate | K4 (k = 3–10) | dimer (k = 3–12) | P4 (k = 10–40) |
|---|---|---|---|
| Σ₁ | 12.3 / 13.1 / 14.1 | 17.2 / 17.2 / 20.3 | 4.5 / 4.6 / 4.9 |
| N1 | 40.6 / 40.5 / 38.5 | 40.0 / 39.9 / 38.1 | 5.8 / 5.8 / 5.7 |
| N2 | 56.9 / 56.9 / 58.2 | 58.8 / 59.0 / 60.4 | 29.6 / 29.6 / 29.5 |
| N2′ path | 56.9 / 56.9 / 58.2 | 58.8 / 59.0 / 60.4 | 28.8 / 28.8 / 28.7 |
| N3 | 71.6 / 72.0 / 76.6 | 83.4 / 83.7 / 86.8 | 74.2 / 74.4 / 75.7 |
| **N3x** | **97.8 / 97.9 / 98.5** | **97.9 / 97.9 / 98.7** | **99.0 / 99.0 / 99.3** |

**Checks (verification order of `plan.md`):**
- **Violations: 0.** No certificate falls below an exact block mass in any block of any
  row, and Σ ≥ S holds for every certificate and row. That is 14 rows × 13 certificates.
- **ρ < 1 certified for every walk matrix and every segment.** The certificate is the M-matrix
  pivot test in `mlu_solve`.
- **Reproduction of the float probe of 2026-10-08:**
  - N1, N2 and N3∧N2 equal the probe re-run at full digits to 6 decimals on all 14 rows. The
    tolerance was 0.05.
  - They agree with the printed one-decimal table to within rounding.
- **N3x against the float N3x (1200 points, the probe's diagnostic):** the rigorous sum at
  M = 2000 sits above the float value by +0.003 to +0.010 on the 12 probe rows (tolerance
  +0.02). The gap between M = 1000 and M = 2000 is 0.003–0.014.
- **Reduction:** the bridge-endpoint reduction equals the full system to the last digit on
  four instances (`results/crosscheck.csv`). One of them, a chain of singletons, is not a
  matching.
- **Precision:** a 256-bit rerun of K4 k=10, dimer k=12 and P4 k=40 (closed forms, and N3x at
  M = 1000) changes no log10 value at 6 decimals over 30 values (`results/aggregate_prec256.csv`).
- **Cost:**
  - pass 2 (final): 298 tasks in 1128 s on 6 workers;
  - the P4 k=40 N3x at M = 2000 takes about 120 s per contour;
  - the closed forms take under 3 min for all rows.

## Hypothesis check (`plan.md`)

1. **Reproduction: supported.** The Arb values with rigorous inputs equal the probe's to
   6 decimals. The probe's float inputs were not an artefact.
2. **Rigorous N3x within +0.02 of the float N3x: supported in pass 2, refuted in pass 1.**
   - **Pass 1** enclosed the block resolvent on each segment as |R(z₀)ᵤᵥ| + r/d². That
     absolute pad swamps the door-to-door entries, which are small by cancellation (about 1/d⁴
     across a distant P4 module). The certificate stayed valid but was 0.05–0.13 orders loose
     on K4 and the dimer, and 2.1, 9.1 and 33.8 orders loose on P4 k = 10, 20 and 40
     (`results/aggregate_pass1.csv`).
   - **Pass 2** uses the entrywise Taylor enclosure Σ_{m<4} rᵐ|(R(z₀)^{m+1})ᵤᵥ|, plus a norm
     tail (r/d)⁴/d/(1−r/d). That keeps the door factor.
3. **Rate shares reproduce `opportunities.md`: supported.** The endpoint secants are
   12.3/17.2/4.5 (Σ₁), 40.6/40.0/5.8 (N1), 56.9/58.8/29.6 (N2) and 71.6/83.4/74.2 (N3∧N2).
   Every value is within 1 point of the expected 12/17/5, 41/40/6, 57/59/30 and 72/83/74.
4. **N3x captures a rate share of at least 0.95 on all three families: supported** (0.978,
   0.979 and 0.990 by endpoint secant).

**MR3 is closed before the 2026-10-15 stop rule.**
- N3x is now a rigorous certificate on all 14 rows. It sits 0.6–5.5 orders above S, and at
  most 2.0 orders above the exact-mass floor.
- For comparison, Σ₁ is 5.6–259.5 orders above S, and the best closed form (N3) is 3.3–80.7.
- What N3x is not: an a priori formula. Phase 3 has to state it as a computer-assisted
  certificate from block data, or else derive the analytic small-block resolvent bounds that
  MR3 lists as the alternative. That alternative was not attempted here.

**Reason 2 of the referee, as measured here.**
- In value, only N3x answers it (≤ 5.5 orders on every row).
- In rate, N3 already captures a share of 0.72–0.83 and N3x 0.98–0.99, against 0.05–0.17
  for Σ₁.

## Figures

- `figures/2026-10-08_neumann_14rows_bounds_vs_k.pdf`, **for the manuscript** (with PNG and
  `figure_manifest.json`).
  - Three panels, K4, dimer and P4, each plotting log10 of S, the floor, Σ₁, N1, N2, N2′ path,
    N3 and N3x against k.
  - **Legend labels set in Phase 3 S9 (2026-10-13)** to the manuscript's names (`LABELS` at the
    top of `code/make_figure.py`): N1 → bridge bound, N2 → walk-sum bound, N2′ (path) → chain
    bound, N3 → hybrid bound, N3x → contour certificate, Σ₁ → Combes–Thomas bound. The figure
    was regenerated; the data did not change. It is Figure `F:linhas` of the manuscript.
    figure-reviewer (Phase 5.6) checks that legend and text agree.
  - The width is 6.0 in. The IOP class is not installed here, so the J. Phys. A text width was
    not measured.

## Pointers

- Aggregate data: `results/aggregate.csv` (long format, one row per row × certificate),
  `results/blocks.csv` (per-block bounds against exact masses) and `results/rate_shares.csv`.
- Raw certified upper endpoints (exact rationals): `results/raw_final.json`,
  `raw_final_pass1.json` and `raw_prec256.json`.
- Oracle: `results/probe_oracle.csv`. It is used for reproduction only.
- Full log: `experiment_log.json`, 92 runs, 0 failures.
  - r0001–r0002: smoke.
  - r0003–r0044: pass 1, superseded (`neumann.py` sha256:a2da7668…).
  - r0045–r0086: pass 2, final (sha256:550a3e56…).
  - r0087–r0092: the 256-bit check.
- Code:
  - `code/neumann.py` (certificates), `code/run.py` (pool, log) and `code/aggregate.py` (sums,
    checks, oracle);
  - `code/analysis.py` (rate shares), `code/crosscheck.py` and `code/make_figure.py`.
- Reproduce:
  ```bash
  cd experiments/2026-10-08_neumann_14rows/code
  ../../../.venv/bin/python -I crosscheck.py --M 300
  ../../../.venv/bin/python -I run.py --tag final                    # ~19 min, 6 workers
  ../../../.venv/bin/python -I aggregate.py --tag final --oracle
  ../../../.venv/bin/python -I analysis.py && ../../../.venv/bin/python -I make_figure.py
  ```
- Environment: `env/requirements.txt`, `env/system_info.json`, `env/git_commit.txt` and
  `env/git_status.txt`. The working tree was not clean at start: the untracked toy,
  `submission (2)/` and two PDFs, as `plan.md` anticipated.

## Notes

- **Deviations from `plan.md`, all made before the final pass:**
  1. **Solves.** The toy's dense `arb_mat.solve` plus semipositivity check was replaced by
     `mlu_solve`, a sparse LU without pivoting.
     - Pivots certified > 0 make I − T a nonsingular M-matrix, which is equivalent to ρ(T) < 1.
     - Only non-negative terms are added, so components near 10⁻²⁷⁰ keep their relative
       accuracy. A dense Arb solve bounds the error normwise and lost such components (a
       9·10⁻⁴ relative error at 10⁻⁴⁰ in the first crosscheck).
  2. **Segment enclosure.** The ball covering a whole segment (toy) inflated 5× in the 4×4
     inverse on P4. The resolvent-identity pad of pass 1 was then too loose, as described above.
     The final enclosure is the Taylor one.
  3. **Bridge rows.** The first version overwrote a |V|U row when a vertex carries two bridges.
     That case cannot occur on the 14 rows, where `rig_setup` asserts a matching. It was caught
     by the Stark smoke test, fixed before pass 2, and added to the crosscheck.
  4. **Reproduction oracle.** The probe was re-run at full digits, by importing it read-only,
     instead of only parsing its printed table. A one-decimal table cannot resolve a 0.05
     tolerance.
  5. **Code layout.** `aggregate.py` was added as a separate file.
- On K4 and the dimer, N2′ path equals N2 to display precision. On P4 it is 0.6–2.5 orders
  looser. Pure N3 equals N3∧N2 on every row.
- **Not run:** the analytic small-block resolvent bounds (MR3's alternative route), and
  general graphs of communities. Both are out of the safe scope.
