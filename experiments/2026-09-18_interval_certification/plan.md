# Plan — 2026-09-18_interval_certification

**Date:** 2026-09-18 · **Scale:** canonical · **Mode:** symbolic (rigorous arithmetic)
**Parent:** opportunities.md #1 · WORKFLOW Tarefa 4 · **Delegated finding:** R-C4
**Supersedes:** `plan_mini.md` (archived as `plan_mini.md.archived`, approved 2026-09-18)

## Question
Does every number printed in Tables `Tab:cliques`, `Tab:dimero` and `Tab:p4`, and
every derived parameter displayed in the text of Examples 7.1–7.3, survive rigorous
evaluation with named, reproducible, outward-rounded arithmetic? Is each `κ_tab` the
minimizer of `Σ1` on the decimal grid of step 1e-4?

## Hypothesis
Supported: for all 14 rows, the printed `Σ1`, `F_esp` and `B²` are at least the upper
endpoint of the Arb ball evaluated at the printed `κ_tab`, and each printed value is
the *tightest* 4-significant-figure upward rounding of that endpoint.
Refuted: some printed value lies below the certified upper endpoint, so the table
under-reports a bound. That is a real inconsistency and goes back to Phase 5.5.
The κ grid argmin is secondary. A disagreement there is an erratum, not a validity
failure, because any κ in (0,1) gives a valid bound.

## Variables
- Independent: the family (K4, dimer, P4) and each printed row `(k, κ_tab)`. For the
  κ audit, all 9999 grid points `κ = j·10⁻⁴`.
- Dependent: balls for `δ, g, W_i, η, ĝ, ϱ_i, γ_κ, c_κ, A_i, L_i, Σ1, F_esp, b_η, B_E,
  B, B²`. The derived columns are `ok_upper`, `ok_tight` and `kappa_is_grid_argmin`.
- Controls: s and t as in §7, and 256-bit working precision (a 512-bit rerun confirms
  the verdicts do not depend on the precision).

## Method
- **Library:** python-flint 0.9.0 (Arb/FLINT 3 ball arithmetic). Each operation
  returns a midpoint–radius ball that rigorously contains the exact result, which is
  outward rounding in ball form. The §7 sentence changes accordingly from "interval
  arithmetic with directed rounding" to "ball arithmetic (Arb via python-flint 0.9.0,
  256-bit precision)". That wording goes to the next address-review wave.
- **Inputs are exact.** `w, Δ, J` are rationals entered as `arb(fmpq)`, never as
  floats. `√26`, `√5`, `φ` and `π` come from Arb. `W_i` and `g` come from the closed
  forms in §7. Their block-spectrum premise (e.g. `spec(P4) = {±φ, ±φ⁻¹}`) is checked
  separately in exact arithmetic with python-flint `fmpz_poly` (characteristic
  polynomial factorization), so no floating eigen-solver enters the certified chain.
- **Distances** are integers from BFS on the explicit graph and are also compared
  with the closed forms `2i`, `2(k−1−i)`, `4i`, and so on.
- **Formulas:** `η = δ/(g−δ)` (eq:eta), `ĝ = (g−2δ)/(2dμ)`,
  `ϱ_i = 1 + 2(W_i+2δ)/(π(g−2δ))`, `γ_κ = log(1+κĝ)`, `c_κ = 1/(1−κ)`, `A_i`, `L_i`
  (eq:Li), `Σ1` (eq:Ti), `F_esp = min{1, Σ1²}`, `b_η` (Cor. comparacao),
  `B_E² = 1 − (h_s−h_t)²/(spread(h)+2dμ)²` (Prop. energia) and
  `B = min{1, Σ1, b_η, B_E}` (eq:Bcert).
- **`min` of balls:** when two balls overlap, the upper endpoint of the min is
  `min(upper_a, upper_b)`, which is still rigorous for an upper-bound claim. Overlaps
  are logged.
- **Tightness:** `printed == ceil₄(upper)`, where `ceil₄` rounds up to 4 significant
  figures in exact rational arithmetic. A printed value that is a valid upper bound
  but not the tightest one is logged as cosmetic (`ok_upper` true, `ok_tight` false).
- **Displayed "≈" parameters** (`η ≈ 0.0050251`, `ϱ_i ≈ 1.135041`, ...) are
  approximations, not bounds. Check that the ball rounds to the printed digits, to
  nearest.
- **κ audit:** evaluate a `Σ1` ball at every grid point, and report the grid argmin
  set (every κ whose lower endpoint does not exceed the minimum upper endpoint).
  `κ_tab` passes if it belongs to that set. When `Σ1` is flat on a plateau (every
  `L_i` taken at `b_i`), 0.5 passes as documented in §7.

## Methodological critique of the plan-mini
1. **The κ objective is undocumented.** §7 says "search on the decimal grid" but not
   of what. `Σ1` is the only printed column that depends on κ (`B` inherits the
   dependence only through `Σ1`), so the objective is taken to be `Σ1`. This
   assumption is logged. If the argmin disagrees, the search is repeated on `B²` before
   anything is flagged.
2. **Certifying `W_i, g` through closed forms alone would be circular** if the closed
   form were wrong. The exact characteristic-polynomial check of each block closes
   that gap.
3. **Not in scope:** Prop. 8.1 and the `C_1024 □ K_64` numbers (Phase 5.5 round 3
   already covered those). Only the tables and the parameters of §7 are in scope.

## Code reuse
No reuse of project code. `verification/build_from_scratch.py` and
`reviews/.../verificar_calculos.py` (Decimal at 70 digits) are read-only
comparators. After the certified run, `code/crosscheck.py` runs them as subprocesses,
never imports them, and records agreement to 1e-12 relative in the log.

## New functions (`code/`)
- `families.py`: `build(family, k) -> (n, edges, h, parts, s, t)`, plus the exact
  rational parameters of each family.
- `certify.py`: `params_ball(fam)`, `sigma1_ball(fam, k, kappa)`,
  `b_ball(fam, k, kappa)`, `ceil4(x: arb) -> fmpq`, `check_row(...)`.
- `block_spectra.py`: exact char-poly factorization of each module (`K4`, the dimer
  block, `P4`).
- `kappa_audit.py`: grid argmin with ball-aware ties.
- `run.py`: parses the three tables straight from `manuscript/main.tex`, read-only,
  so the check always targets the current text. Writes `results/aggregate.csv`,
  `results/params.csv` and `results/kappa_audit.csv`, and appends to
  `experiment_log.json`.
- `crosscheck.py`: runs the comparators as described above.

## Verification order
1. Exact block spectra match §7 (`spec(K4 + Δi)`, `24i+5±√26`, `±φ, ±φ⁻¹`).
2. BFS distances match the closed forms, and `D` matches the table column.
3. (H1): `‖H_out‖ = w`, since the bridges form a matching. Check the matching
   property explicitly. (H2): `g > 2δ`.
4. The displayed parameters round correctly.
5. `Σ1, F_esp, B²` balls against the printed values, giving `ok_upper` and `ok_tight`.
6. The κ audit.
7. The 512-bit rerun gives identical verdicts.
8. Cross-check against the two comparators.

## Files
- **To create:** `code/*.py`, `results/{aggregate.csv, params.csv, kappa_audit.csv,
  viz_schema.json}`, `experiment_log.json`, `README.md`,
  `env/{requirements.txt, system_info.json}`.
- **Reuse read-only:** `manuscript/main.tex`, `verification/`, `reviews/.../verificar_calculos.py`.
- **Update:** none outside this folder. The edits to `main.tex` go to the next
  address-review wave.

## Time
Build ≈ 1 h. The run is 14 rows × 9999 κ balls ≈ 1.4·10⁵ evaluations, a few minutes.
