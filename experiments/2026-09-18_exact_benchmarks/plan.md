# Plan — 2026-09-18_exact_benchmarks

**Date:** 2026-09-18 · **Scale:** canonical · **Mode:** numerical
**Parent:** opportunities.md #1 · WORKFLOW Tarefa 4 · **Delegated finding:** R-B5
**Supersedes:** `plan_mini.md` (approved 2026-09-18)

## Question
How far do the certificates sit from the quantities they bound, on every row of
Tables 3–5? Does the certified spatial decay track the true decay of transport?

## Hypothesis
Supported: `Σ1/S` stays within a few orders of magnitude and changes smoothly with
`k`. The true decay rate `γ_true` (the slope of `−log S` against `D`) is at least the
certified `γ_κ`, and the gap is quantified. The masses `‖P_i e_x‖` follow the profile
of `L_i(x)`.
Refuted: a ratio that grows like `e^{(γ_true − γ_κ)D}` with a large exponent gap. That
does not make the paper wrong, but the conservativeness gets reported openly,
together with its rate.

## Variables
- **Benchmark (all 14 rows):** K4 k ∈ {3,4,5,6,8,10}, dimer k ∈ {3,5,7,9,12},
  P4 k ∈ {10,20,40}, with s and t as in §7 and `κ = κ_tab`.
- **Figure sweep (K4 family, w = 0.1, m = 4):**
  (a) k = 3…12 at Δ = 24, versus `D = 2k−1`;
  (b) k = 5 fixed, with Δ chosen so that `g/(2δ) = (Δ−4)/0.2` takes the values
  {100, 30, 10, 5, 3, 2, 1.5, 1.2, 1.05}, down to the edge `g → 2δ⁺`. In (b), κ is
  optimized on the grid in double precision. The figure is illustrative, not a table.
- **Dependent:** exact `S`, ITF `S²`, exact `F̄ = Σ c_λ²`, and the `F_max` bracket
  `[F_samp, S²]`. The masses are `‖P_i e_s‖` and `‖P_i e_t‖` for every block i,
  compared with `L_i`. The ratios are `Σ1/S`, `B²/S²` and `F_esp/F_samp`, and
  `γ_true` comes from a least-squares fit per family.

## Method
- **Precision:** mpmath `eigsy` at 60 digits for every benchmark row. Measured cost:
  n=40 takes 1.7 s and n=80 takes 12.5 s, so n=160 should take about 100 s. The
  certificates go down to 1e-22, and double-precision eigenvectors cannot resolve
  projector entries below about 1e-16·‖H‖, so double precision serves only as a check
  on rows where S > 1e-10. A 40- versus 60-digit rerun on the P4 k=40 row confirms
  convergence to 3 significant figures.
- **Distinct eigenvalues:** eigenvalues are grouped with tolerance 1e-40, and the log
  records the smallest gap between groups, which must be ≫ 1e-40. The exact
  degeneracies of K4 (the unbridged vertices keep an S2/S3 symmetry) therefore merge
  correctly.
- **Block projectors:** `P_i = Σ_{λ ∈ I_i} Π_λ`, with `I_i = [min Σ_i − δ, max Σ_i + δ]`.
  Every eigenvalue must land in exactly one `I_i`, with `|V_i|` eigenvalues each,
  which is a direct check of Lemma aglo.
- **F_samp:** the supremum over τ ∈ [0, T] with T = 2π·10³/λ_min-gap, sampled at a
  step of 0.05/‖H‖, followed by golden-section refinement at the top 20 local maxima.
  It is a lower bound for `F_max`. The sup over τ ≥ 0 is not attainable, and the
  bracket says so.
- **Separation of use cases:** the certificates come from the Arb pipeline of
  `2026-09-18_interval_certification` (read from its `results/aggregate.csv`), never
  from the diagonalization. The README states that the exact values exist only for
  calibration.

## Methodological critique of the plan-mini
1. **k=3 is structurally easy:** there, `B²` is set by `b_η²` and not by `Σ1`. It
   stays as the entry point the card asks for, and the most informative rows are the
   ones where `Σ1` is active (K4 k≥5, dimer k≥5, all of P4). The README says so.
2. **`F_max` cannot be computed exactly.** With a quasi-periodic amplitude, sampling
   gives only a lower bound. The bracket `[F_samp, S²]` is the honest object, and
   `S²` is the ITF, which is the quantity the certificate actually bounds. So the
   ratio headline is `Σ1²/S²`, and `F_samp` is supplementary.
3. **Confounder in (b):** as Δ drops, the site-energy range shrinks, so `B_E` becomes
   competitive and `b_η` grows. Panel (b) plots all bounds separately and makes no
   claim that `Σ1` alone is best.

## Code reuse
None from the project. The primitives are written fresh in `code/spectral_mp.py`
(mpmath) and `code/families.py`, which is copied from the certification experiment
with its SHA-256 logged. `verification/` is not imported.

## New functions
- `spectral_mp.py`: `eig_groups(H, dps, tol) -> (lams, Pis)`,
  `block_projectors(lams, Pis, intervals)`, `S_exact`, `Fbar_exact`,
  `F_sampled(c, lams, T, step)`.
- `bounds.py`: double-precision `Σ1, b_η, B_E` for the sweep only, cross-checked
  against the Arb values on the table rows.
- `run_benchmark.py`, `run_sweep.py`, `make_figure.py` (SciencePlots and the
  figures_style rcParams, at LAA single-column width).

## Verification order
1. Lemma aglo counts per `I_i` are correct on every instance.
2. `Σ_i P_i = I` and `P_i² = P_i` to 1e-50.
3. `S ≤ Σ1`, `F̄ ≤ Σ2`, `F_samp ≤ S²` and `‖P_i e_x‖ ≤ L_i(x)` on every row. Any
   violation is a real inconsistency and stops the run.
4. The 40- versus 60-digit convergence.
5. Only then compute the ratios and fit `γ_true`.

## Files
- **To create:** `code/*`, `results/{aggregate.csv, masses.csv, sweep.csv,
  viz_schema.json}`, `figures/{itf_vs_bounds.pdf, figure_manifest.json}`,
  `experiment_log.json`, `README.md`, `env/*`.
- **Reuse read-only:** `../2026-09-18_interval_certification/results/aggregate.csv`.

## Time
Build ≈ 1.5 h. Run: the benchmark takes about 10 min (P4 k=40 dominates), and the
sweep takes a few minutes.
