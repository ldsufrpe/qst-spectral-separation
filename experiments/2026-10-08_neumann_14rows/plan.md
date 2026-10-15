# Plan — 2026-10-08_neumann_14rows

**Date:** 2026-10-08 · **Scale:** canonical · **Modes:** numerical (benchmark against exact values)
and symbolic (certified interval computation, Rules Sym3/Sym12)
**Parent:** opportunities.md #2 (chosen), MR2 and MR3 · `derivacoes/2026-10-08_neumann_pontes/`
**Calibrated by:** `experiments/toy/2026-10-08_neumann_calibration/` (T1, T2). That toy's
plan-mini and plan were approved by the author on 2026-10-08, and they stand in for this
experiment's Stage 1A.
**Target venue:** J. Phys. A: Mathematical and Theoretical.
**MR3 stop rule:** 2026-10-15. If the rigorous N3x is not done by then, it stays a diagnostic,
and the README says so.

## Question
On the 14 rows of Table 3, how far do the Neumann-route certificates sit from the exact S? The
certificates are N1, N2, N2′, N3 and a **rigorous** N3x. Distance is measured two ways: in value,
log10(Σ/S), and in rate, as the share of the exact decay rate in k (MR2). The question also covers
whether N3x, computed rigorously, keeps the probe's distance of 10^0.6–10^3.2, which would close
MR3.

## Hypothesis
- **Supported:**
  1. The Arb evaluation with rigorous inputs reproduces the probe's log10 Σ (`sondagem_n123_saida.txt`)
     to within ±0.05 on every row and bound. The toy reproduced it to 10⁻⁶ on 3 rows.
  2. The rigorous N3x upper sum at M = 2000 sits within +0.02 in log10 of the float N3x of
     `sondagem_portas_saida.txt` on the 12 rows that the probe covered. It also produces the
     first values for P4 k=20 and k=40.
  3. The rate shares of `opportunities.md` reproduce under the endpoint-secant definition, to
     within ±3 percentage points:

     | bound | K4 | dimer | P4 |
     |---|---|---|---|
     | Σ₁ | 12% | 17% | 5% |
     | N1 | 41% | 40% | 6% |
     | N2 | 57% | 59% | 30% |
     | N3 | 72% | 83% | 74% |
  4. N3x captures ≥ 95% of the exact rate on all three families.
- **Refuted, and what each refutation means:**
  - A row that moves by more than 0.05 means the probe had a float artefact.
  - Rigorous N3x more than 0.1 above the float value means the quadrature was hiding loss.
  - An N3x rate share below 90% on P4 means MR3 answers reason 2 on K4 and the dimer only.
  - **Any bound below an exact mass, or any Σ < S, is a correctness failure.** It stops the
    run and goes back to Phase 5.5 material.

## Variables
- **Instances:** K4 k ∈ {3,4,5,6,8,10}, dimer k ∈ {3,5,7,9,12}, P4 k ∈ {10,20,40}. s and t are as
  in §7, and the families are built exactly by `families.build`.
- **Certificates.** For each one, the per-block Lᵢ(x) for x ∈ {s,t}, and Σ = Σᵢ Lᵢ(s)Lᵢ(t):
  - Σ₁, the paper's Combes–Thomas certificate. It is read from
    `2026-09-18_interval_certification/results/aggregate.csv` and is not recomputed.
  - N1 (Theorem N1).
  - N2 (Prop. N2, the community-graph walk sum).
  - N2′ in two forms:
    - as stated (sorted distances G₍ₘ₎);
    - the ladder form q^N/(2N−1)!!, only where g_il ≥ |i−l|g is certified.

    On these rows every bridge has the same weight w = δ, so the per-bridge variant coincides with
    the global one.
  - N3 (Prop. N3), minimized over κ ∈ {1/2, 7/10, 4/5, 9/10, 19/20}. Minimizing over κ is
    legitimate because the bound holds for each κ. The probe's N3 was min(N3, N2), and that
    combination is reported separately.
  - N3x as a rigorous upper Riemann sum at M ∈ {1000, 2000} segments per contour. The minimum of
    the two is also a valid certificate.
  - "best": the per-block minimum over all of the above, a valid certificate in its own right.
- **References,** from `2026-09-18_exact_benchmarks`:
  - exact S;
  - exact masses ‖Pᵢeₓ‖;
  - floor = Σᵢ‖Pᵢeₛ‖‖Pᵢeₜ‖.
- **Dependent:**
  - log10(Σ/S) and log10(Σ/floor);
  - the rate share per family, computed three ways:
    - **endpoint secant,** first and last k (the `opportunities.md` definition);
    - **least-squares slope** of log10 Σ against k over all rows;
    - **last-pair secant** (the asymptotic proxy).

    The share is slope(Σ)/slope(S). The η floor flattens small k, so the three definitions
    differ, and all three are reported.
  - The N3x quadrature gap, log10 Σ(M=1000) − log10 Σ(M=2000).
  - Wallclock per row.

## Method
- **Rigorous inputs,** as in the toy:
  - exact characteristic polynomials and Arb root enclosures for the block hulls;
  - exact rational safe-side g_ij, g (lower), Wᵢ (upper) and δ = w, exact because the bridges form
    a matching (asserted);
  - η = δ/(g⁻−δ).

  Every walk matrix T has ρ(T) < 1 certified by semipositivity: a point x > 0 with (I−T)x > 0,
  verified in Arb.
- **N3x, rigorous.** The stadium sits around the rational outer enclosure of Jᵢ with ρ = g⁻/2,
  and the Riesz conditions are checked exactly. On each segment:
  - the block resolvents are enclosed with `acb_mat.inv` on a ball that covers the segment;
  - U is the entrywise upper endpoints of the moduli;
  - the bound is ‖Qᵢ U(I−|V|U)⁻¹eₓ‖ times the upper length.

  **Reduced system.** Only the bridge endpoints B, with |B| = 2(k−1), carry unknowns:
  (I − |V|_BB U_BB) y_B = (eₓ)_B + [x∉B]|V|_BB U_Bx, and then (Qᵢm)_u = Σ_{v∈B∩Vᵢ} U_uv y_v. The
  spectrum of |V|U equals that of |V|_BB U_BB plus zeros, so certifying ρ < 1 on the reduced
  matrix suffices.
- **Precision:** 128 bits for the run. Three rows (K4 k=10, dimer k=12, P4 k=40) are rerun at
  256 bits for N1/N2/N3, and their N3x is rerun at M=1000.
- **Parallelism:** a `multiprocessing` pool of 6 workers over (row, block i) tasks. Arb is
  deterministic, so the results do not depend on scheduling. `OPENBLAS_NUM_THREADS=1`.
- **Rate shares:** computed by `analysis.py` from `aggregate.csv` only.

## Methodological critique (Stage 1B)
1. **The rate share is not the same as rate sharpness.** On P4 there are only 3 points (k = 10,
   20, 40), and S decays superexponentially, so the share depends on the window. That is why
   three definitions are reported and none is singled out. The headline in the README is the
   endpoint secant, for continuity with `opportunities.md`, and the last-pair secant is shown
   next to it.
2. **N3x is a per-instance computer-assisted certificate, not an a priori formula.** It uses
   only block data: block resolvents along Γᵢ, δ and the bridges. It is rigorous, but it is not a
   closed-form bound like N1–N3. The README must say this plainly, because Phase 3 has to decide
   how to state it. The "analytic bounds on small-block resolvents" route of MR3 is not
   attempted. It would give a formula, and it stays as a Phase 3 option.
3. **κ minimization and the "best" column can look like cherry-picking.** Each is a valid bound
   for every choice, so a minimum over choices fixed in advance is itself a bound. The κ grid is
   fixed here, before the run.
4. **Σ₁ is not recomputed.** It is read from its own certified experiment. The comparison is
   fair only because both sides are Arb upper bounds on the same S with the same s and t, and
   that is asserted.

## Code reuse (read-only imports, never edited)
- `2026-09-18_exact_benchmarks/code/families.py`: `build`, `bfs`.
- `2026-09-18_interval_certification/code/certify.py`: `block_spectrum` (l.74), `upper_q`,
  `lower_q`, `amin`.
- `2026-09-18_exact_benchmarks/results/{final,masses}.csv`: S, masses.
- `2026-09-18_interval_certification/results/aggregate.csv`: Σ₁.
- `derivacoes/2026-10-08_neumann_pontes/sondagem_{n123,portas}_saida.txt`: the probe's values,
  used only as the reproduction oracle. They are parsed, and the scripts are not re-imported.
- The toy's `code/run.py`: `rig_setup`, `certify_M`, `bounds_123`, `stadium_segments` and
  `n3x_upper` are copied into `code/neumann.py`, with the toy file's SHA-256 logged.
- `2026-09-18_exact_benchmarks/code/logutil.py`: copied, with the header updated to
  opportunity #2 and J. Phys. A, and its SHA-256 logged.

## New code
- `code/neumann.py`:
  - the toy functions, plus `n2prime(P, form)`;
  - `n3x_upper_reduced(P, i, M)`;
  - `ladder_certified(P)`.
- `code/run.py`: the row loop with the pool, appending to `experiment_log.json`. It runs the
  `--smoke` rows (K4 k=3, M=200) or `--tag final`.
- `code/crosscheck.py`: the reduced N3x against the dense N3x of the toy on K4 k=3 and dimer k=3,
  at the same M. They must agree to 10⁻²⁵ relative.
- `code/analysis.py`: rate shares and the three slope definitions, written to
  `results/rate_shares.csv`.
- `code/make_figure.py`: `figures/neumann_vs_exact.pdf`, three panels (K4, dimer, P4) of log10 Σ
  against k for every certificate, with S and the floor. Style follows `figures_style.md` at
  J. Phys. A width, and `figure_manifest.json` is written alongside.

## Verification order
1. **Inputs:** (H2) holds rigorously on all 14 rows, and the ladder certification is logged
   per row.
2. **Smoke:** K4 k=3, all bounds, M = 200, under 30 s.
3. **Reduction:** `crosscheck.py` must pass, reduced = dense, before any large-row N3x.
4. **Reproduction:** N1/N2/N3∧N2 against the probe, ±0.05, on all 14 rows. Any miss is
   investigated before going on.
5. **Validity:** every bound is ≥ its exact mass, and Σ ≥ S, for every certificate and row. ρ < 1
   is certified everywhere.
6. **N3x:** the upper sum is ≥ the float N3x on the 12 probe rows, and M=2000 ≤ M=1000. The second
   is expected rather than guaranteed, and any breach is logged.
7. **Precision:** the 256-bit rerun agrees.
8. **Only then:** metrics, rate shares and the figure.

## Files
- **To create:**
  - `code/{neumann.py, run.py, crosscheck.py, analysis.py, make_figure.py, logutil.py}`;
  - `results/{aggregate.csv, blocks.csv, rate_shares.csv, crosscheck.csv, viz_schema.json}`;
  - `figures/{neumann_vs_exact.pdf, figure_manifest.json}`;
  - `env/{requirements.txt, system_info.json, git_commit.txt, git_status.txt}`;
  - `experiment_log.json`, `README.md`, `report_validation.md`.
- **To update:** none outside this folder. WORKFLOW.md gets a status line at the end of the
  session, as is the project's convention.
- **Read only:** everything under Code reuse. `manuscript/main.tex` is not touched.

## Stopping criterion
All 14 rows × all certificates, plus N3x at M ∈ {1000, 2000}. If one row's N3x at M=2000 exceeds
2 h of wallclock, keep M=1000 for that row and log the reason.

## Acceptance preview (J. Phys. A)
- A certified table:
  - every value is an Arb upper bound, checked against exact values computed independently;
  - zero violations;
  - reproducible from one command.
- MR2 is satisfied by the rate-share table.
- MR3 is closed if N3x is rigorous on all 14 rows, with the honest label "computer-assisted
  certificate from block data".

## Anticipated risks
- **N3x cost on P4 k=40.** With the reduced dense solve of size 78, the estimate is about
  30 ms per segment, so 40 contours × 2000 segments ≈ 40 min on 6 workers. Mitigation: the M=1000
  fallback.
- **Enclosure blow-up near the P4 contours.** ρ = g/2 ≈ 0.88 and the segments are short, so this is
  not expected. The toy found none on K4 or the dimer.
- **Git state.** The working tree has untracked files: the toy, `submission (2)/`, and two PDFs.
  This is recorded in `env/git_status.txt`. The experiment code is self-contained, and the toy is
  committed with it at the end.

## Time
- Build: about 2 h.
- Run: smoke < 1 min; N1–N3 on all rows < 2 min; N3x about 1 h wall on 6 workers; analysis and
  figure, minutes.
- Stage 4 (README, prose-author, gate): about 1 h.
