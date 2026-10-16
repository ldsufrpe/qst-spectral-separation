# Plan — 2026-10-15_rc3_threshold

**Date:** 2026-10-15 · **Scale:** canonical · **Modes:** numerical (exact S in multiprecision)
and symbolic (certificates in Arb, verdicts in exact rationals)
**Parent:** card **R-C3** of the Phase 6 review of opportunity #2 (`peer_review_reconciled.json`;
`.sciflow/review_triage.json`, decision "accepted", delegated to experiment-runner on 2026-10-15).
This experiment closes that delegation.
**Target venue:** J. Phys. A: Mathematical and Theoretical.
**Stage 1A:** `plan_mini.md.archived`, approved by the author on 2026-10-15 with the recommended
choice on all four decisions. Those choices are the stable threshold k\*, the extension on K4 and
the dimer only, six values of q, and no figure by default.
**Pre-plan probe:** `exploration/preplan_probe/`. It used M = 300, κ step 1/50 and float64 S. It
informed the design, and none of its numbers is reported.

## Question
For each family of `Tab:familias`, at which chain length k does each certificate guarantee
Σ₁ < √2 − 1, hence F_Bloch < 2/3 at every time (`cor:fent`)? The certificates are the bridge
bound, the walk-sum bound, the chain bound, the contour certificate and the Combes–Thomas
baseline. The same question is asked of the exact overlap sum S. On K4 and the dimer every
certificate crosses already at k = 2, so the question is repeated there on stronger bridges.

## Hypothesis
- **Supported:**
  - **P4 at the published parameters:** k\* = 2 for the contour, walk-sum and chain bounds and for
    S; k\* = 4 for the bridge bound; k\* = 9 for the baseline. The refined certificates then
    change the baseline's verdict on k = 3…8, at μ = 1.
  - **K4 and the dimer:** some q of the grid gives a k ≤ 10 at which the baseline fails and the
    contour certifies. For q ≥ 0.35, 2η > √2 − 1, so at k = 2 only the contour certifies, as S
    does.
- **Refuted:**
  - no q of the grid gives such a k on K4 or the dimer. Then R-C3 takes its second route, and
    the numerical claim is delimited to quantitative suppression;
  - a verdict differs from the probe on a row whose Σ₁ lies more than 1 % from √2 − 1. That is
    investigated before anything else.
- **Correctness failure (stop):**
  - any Σ₁ < S;
  - any per-block bound below its exact mass;
  - a Lemma `aglo` count that fails.

## Variables
- **Published rows** (module, Δ and w as in `Tab:familias`):
  - K4: k = 2, 3;
  - dimer: k = 2, 3;
  - P4: k = 2…10.

  The upper end is the first row of `Tab:linhas`.
- **Stronger bridges,** on K4 and the dimer, k = 2…10. Only the bridge weight w changes. The
  grid is q = 2w/g ∈ {0.20, 0.25, 0.30, 0.35, 0.40, 0.45}:
  - K4 (g = 20): w = 10q, exactly 2, 5/2, 3, 7/2, 4, 9/2;
  - dimer (g = 24 − 2√26): w = q·g/2 rounded to two decimals, 138/100, 173/100, 207/100,
    242/100, 276/100, 311/100. The CSV records the q that each rational w actually gives.

  Every grid point satisfies g > 2δ and q′ = q < 1/2. The bridges form a matching, so δ = w and
  q′ = q.
- **Per row:**
  - log10 Σ₁ for each of the five certificates, from the exact rational upper endpoint;
  - log10 S;
  - the verdict of each one against √2 − 1;
  - η, q and μ.
- **Per (family, w, certificate):**
  - k_first, the first k that certifies;
  - **k\***, the smallest k from which it certifies at every scanned k. On the published
    families, the rows of `Tab:linhas` beyond the scan also have to certify. They are read from
    `2026-10-08_neumann_14rows/results/aggregate.csv`, not recomputed;
  - the list of k that fail.

  When the certificate fails at k = 10, k\* is reported as "> 10".
- **Fixed, as in `Tab:linhas` and Appendix A:**
  - contour with M = 2000 arcs, at 128 bits;
  - CT with κ minimizing Σ₁ on the grid κ = j/10⁴, j = 1…9999;
  - s and t at the two ends of the chain.
- **Onset (decided by rule, before the run):** for K4 and the dimer, the onset is the smallest q
  of the grid with a k ≤ 10 at which CT fails and the contour certifies. The `Threshold` macros
  report the onset and q = 0.40.

## Method
- **Certificates.** Each is Σ₁ of `teo:principal` with s and t at the two ends. Each per-block
  factor is the upper endpoint of an Arb ball, as an exact rational. The factors come from:
  - bridge: Theorem `teo:pontes`, `bounds_123` N1;
  - walk-sum: Lemma `prop:passeios`, `bounds_123` N2;
  - chain: Theorem `cor:cadeia`, `eq:cadeiamassa` with the coupling of each bridge, `n2prime`
    N2p_path;
  - contour: Proposition `prop:certcontorno`(b), `n3x_upper_reduced` at M = 2000;
  - CT: the baseline of §7, `certify.sigma1`, minimized over the κ grid. A minimum over a grid
    fixed in advance is itself a bound.
- **Verdict of a certificate.** The test (Σ₁ + 1)² < 2 is exact in rationals.
- **Exact S.** `spectral_mp` at adaptive precision, with the rule of `run_benchmark.py`:
  - start at 60 digits;
  - raise until S > 10^−(dps−30);
  - if the precision was raised, confirm at dps + 20 to 10⁻⁶ relative.

  The exact block masses, on the intervals I_i = J_i′ + [−δ, δ] of `rig_setup`, feed the
  validity check. The verdict on S is numerical, so its margin to √2 − 1 is reported.

## Methodological critique (Stage 1B)
1. **The published K4 and dimer rows pass for a structural reason.** η is 0.005 and 0.004, and
   Σ₁ ≤ 2η + (k − 2)η² sits far below the threshold whatever the certificate. They are kept as
   the literal answer to R-C3, and the extension carries the content. P4 at k = 2 also passes on
   the η floor alone (2η = 0.4099, 1 % below 0.4142). The stable threshold k\* keeps that from
   reading as a verdict.
2. **At k = 2, every norm-based bound equals 2η.** With N = 1 the factor ϱq/((1−q)√(1−η²)) of the
   bridge, walk-sum and chain bounds exceeds η, so all three are capped. Only the contour, which
   reads the entries of the block resolvents, can go below η. On the extension, where 2η > √2 − 1,
   k = 2 therefore separates the contour from the rest by construction. The README states this,
   so that the reader does not take it for a numerical accident.
3. **The baseline pays twice for a strong bridge.** Lemma `lema:ct` needs μ ≥ every coupling,
   bridges included. On the extension w > 1, so μ = w, and the CT rate per edge
   log(1 + κ(g − 2δ)/(2dw)) also falls with w. This is the lemma as stated, not a choice here, and
   it confounds the extension. The clean evidence is therefore P4 at its published parameters,
   where μ = 1 and the baseline still fails at k = 3…8. The CSV records μ on every row, and the
   README names the confound.
4. **S, not the supremum of the amplitude.** `cor:fent` needs a bound on |f(τ)|, and S ≥ sup|f|.
   The first crossing of S is therefore an upper estimate of the true first crossing. The probe
   puts S ≤ 7·10⁻³ at k = 2 everywhere, so the exact column is 2 and the estimate is exact here.
5. **k\* depends on the scanned range.** It is "stable within k ≤ 10", plus the tabulated rows on
   the published families. Nothing is claimed beyond k = 10 on the extension, and the README says
   so.
6. **Fairness.** All five certificates use the same s, t, δ, g and η, the same exact inputs and
   the same verdict test. CT keeps its κ freedom, as in `Tab:linhas`.

## Code reuse (read-only imports, never edited)
- `2026-10-08_neumann_14rows/code/neumann.py`:
  - `rig_setup` (l.119);
  - `bounds_123` (l.212);
  - `n2prime` (l.295);
  - `n3x_upper_reduced` (l.488);
  - `lg10` (l.42).

  `rig_setup` builds the family through the module-level name `build` (l.30). The new code
  rebinds `neumann.build` inside a scoped context manager, which leaves the file untouched. Each
  pool worker runs its tasks one at a time, so the scope is safe.
- `2026-09-18_interval_certification/code/certify.py`:
  - `parameters` (l.114) and `geometry` (l.127). Both read G["w"], G["d"] and G["mu"], so a
    re-weighted G works unchanged;
  - `sigma1` (l.133);
  - `upper_q` and `lower_q` (l.40, l.45).
- `2026-09-18_exact_benchmarks/code/families.py`: `build` (l.57) and `dense` (l.114). It is
  identical to the copy in `interval_certification`, which `neumann.py` resolves first, as checked
  with `diff`.
- `2026-09-18_exact_benchmarks/code/spectral_mp.py`: `to_mp` (l.15), `eig` (l.24), `residuals`
  (l.37), `group` (l.49), `coefficients` (l.64), `S_Fbar` (l.69) and `block_masses` (l.73).
- `2026-09-18_exact_benchmarks/code/run_benchmark.py`: the adaptive rule (l.102–124) and `mpq`
  (l.31) are **copied**, not imported. Importing that module would pull in its own `logutil` under
  the same module name.
- `2026-10-08_neumann_14rows/code/logutil.py` is copied with the header and provenance updated,
  and its SHA-256 is logged. The pool, checkpoint and logging pattern follows
  `2026-10-08_neumann_14rows/code/run.py`.
- **Reference values** (read only):
  - `2026-10-08_neumann_14rows/results/aggregate.csv`;
  - `2026-09-18_interval_certification/results/aggregate.csv`;
  - `2026-09-18_exact_benchmarks/results/final.csv`.

## New code
- `code/rc3.py`:
  - `build_w(name, k, w)`: `families.build`, then the bridge weights replaced and μ recomputed;
  - `setup(name, k, w)`: `rig_setup` with `build_w`, plus asserts of (H2), q′ < 1/2 and δ = w;
  - `closed(P, G)`: per-block bridge, walk-sum and chain factors, and `ct_best(G)`, which
    returns (Σ₁, the κ argmin set);
  - `contour(P, i, M)`;
  - `exact(G, P, dps0)`: S, the masses, the dps trail and the checks;
  - `certifies(sig)`: the exact test;
  - `thresholds(rows)`: k_first, k\* and the failures.
- `code/run.py`: tasks ("closed", fam, w, k), ("contour", fam, w, k, i) and ("exact", fam, w, k),
  on a pool of 6 workers with `OPENBLAS_NUM_THREADS=1`. It checkpoints `results/raw_<tag>.json`
  after every task and writes one log entry per (family, w, k). It also takes `--smoke`
  (K4 k = 3, published w, M = 200).
- `code/aggregate.py`. It writes:
  - `results/aggregate.csv`, one row per (family, w, k, certificate) with S as a certificate
    "exact";
  - `results/thresholds.csv`, one row per (family, w, certificate);
  - `results/reproduction.csv`;
  - `results/viz_schema.json`.

  It also runs the validity checks.
- **No manuscript figure,** by the 1A decision. `viz_schema.json` is the dashboard contract
  only. If the results ask for a figure, that comes back to the author first.

## `Threshold` macros (scheme fixed now)
`manuscript/make_experimental_macros.py` gains a block that reads `results/thresholds.csv` and
`aggregate.csv`. It writes:
- **P4 at the published parameters:**
  - `\ThresholdPfour{Exact,Contour,Chain,Walk,Bridge,CT}`, the k\* values;
  - `\ThresholdPfourCTFailMin` and `\ThresholdPfourCTFailMax`, the range of k where CT fails and
    the contour certifies;
  - `\ThresholdPfourBridgeFail`;
  - `\ThresholdPfourEtaFloor`, which is 2η.
- **K4 and the dimer at the published parameters:** `\ThresholdPublishedSigmaMax`, the largest Σ₁
  over every certificate and row with k ≤ 3, which shows that everything certifies there.
- **Onset and q = 0.40, for K4 (`Cliques`) and the dimer (`Dimer`):**
  - `\Threshold<Fam><Onset|Strong>{Q,W}`;
  - the six k\* values `\Threshold<Fam><Onset|Strong>{Exact,Contour,Chain,Walk,Bridge,CT}`;
  - `\Threshold<Fam><Onset|Strong>CTFail{Min,Max}`.

  A k\* that is "> 10" is written as `$>10$`. `\ThresholdScanMax` is 10.

Every macro carries a comment naming its CSV and column. The prefix is `Threshold`. The script
asserts that no other macro changes, by comparing the regenerated file outside the new block with
the current one.

## Verification order
1. **Inputs.** `build_w(name, k, w_published)` equals `families.build(name, k)` on edges, h, μ
   and d for all three families. (H2) and q′ < 1/2 hold exactly on every row. η and q match
   their closed forms.
2. **Smoke.** K4 k = 3, published w, all five certificates and S, M = 200. It runs in under a
   minute.
3. **Reproduction,** on the published rows that `Tab:linhas` already has (K4 k = 3, dimer k = 3,
   P4 k = 10):
   - bridge, walk-sum, chain and contour (M = 2000) against the `log10_sigma` of
     `2026-10-08_neumann_14rows/results/aggregate.csv` (N1, N2, N2p_path, N3x_M2000), to 10⁻⁶;
   - CT against `Sigma1_grid_min_upper` of `interval_certification`, to 10⁻⁶ relative;
   - S against `final.csv`, to 10⁻⁶ relative.

   A miss stops the run.
4. **Validity,** on every row:
   - Σ₁ ≥ S for every certificate;
   - every per-block factor ≥ its exact mass;
   - the `aglo` counts, the residual ≤ 10^−(dps−10), and dps convergence where it was raised;
   - ρ < 1 certified for the walk matrices and for every contour arc (`NotMMatrix` fails the task).
5. **Verdicts** from the exact rational test. On every row whose Σ₁ lies more than 1 % from
   √2 − 1, the verdicts agree with the probe.
6. **Thresholds:** k_first, k\* and the failures, and the onset by the rule above.
7. **Macros.** Regenerate them. `main.tex` compiles with exit 0, and the `pdftotext` of
   `main.pdf` is identical before and after, because `main.tex` uses no `Threshold` macro yet.
   `audit_citations` is unchanged.

## Files
- **To create** (under `experiments/2026-10-15_rc3_threshold/`):
  - `plan.md`, with `plan_mini.md.archived` kept;
  - `code/{rc3.py, run.py, aggregate.py, logutil.py}`;
  - `results/{raw_smoke.json, raw_final.json, aggregate.csv, thresholds.csv, reproduction.csv,
    viz_schema.json, run_final.log}`;
  - `env/{requirements.txt, system_info.json, git_commit.txt, git_status.txt}`;
  - `experiment_log.json`, whose `provenance_note` cites R-C3;
  - `README.md`, with the narrative written by the prose-author subagent;
  - `report_validation.md`;
  - `exploration/preplan_probe/{probe.py, probe_output.txt}`, which already exists.
- **To update:**
  - `manuscript/make_experimental_macros.py`, the `Threshold` block;
  - `manuscript/experimental_macros.tex`, regenerated;
  - **the "Pendência da Phase 4", documentation only.** Results and code are untouched, and each
    correction is marked with a dated note:
    - `experiments/2026-09-18_exact_benchmarks/plan.md`, the F_samp item (l.47–50). The code uses
      step π/(8·spread), T = min(200·2π/ω_min, 2·10⁶ steps) and drops terms with |c| < 10⁻⁸ S.
      It refines the 20 best samples with bounded Brent minimization (`spectral_mp.F_sampled`,
      l.83–124);
    - `experiments/2026-10-08_stark_sharpness/README.md` (l.36–37 and l.95). The kept run is the
      more precise one, at 70–214 digits, and it is compared with a run 20 digits below it
      (`code/run.py`, `exact_adaptive`, l.102–114);
  - `WORKFLOW_historico.md`, "R1 feita" at the end;
  - `WORKFLOW.md`: R1 ☑ and the Phase 4 pendency removed.
- **Read only:**
  - everything under Code reuse;
  - `manuscript/main.tex`, which is not edited.

## Stopping criterion
All 13 published rows and 108 extension rows, with the five certificates and S on each, and
zero violations. If a contour row takes more than 30 min, it falls back to M = 1000, and the log
says so. Every value is a valid bound for any M.

## Acceptance preview (J. Phys. A)
- R-C3 gets, per family, the k at which each certificate first certifies F_Bloch < 2/3 for good,
  next to the exact crossing.
- It also gets a regime, published P4 included, where the refined certificate changes the
  baseline's verdict. If the extension fails, the delimitation route applies instead.
- The numbers reach §7 only through the `Threshold*` macros, in session R2.
- Zero violations. The overlapping rows reproduce the published campaigns.

## Anticipated risks
- **Contour cost.** The extrapolation from the 2026-10-08 log (K4 k = 10: 172 s; P4 k = 10:
  170 s at M = 2000) gives ~1.4 h CPU over 121 rows. With CT and the exact S, the total is
  ~2.5 h CPU, which is ~30 min wall on 6 workers.
- **Precision of S** on the extension at k = 10 (S ~ 10⁻²⁵ to 10⁻⁴⁰). The adaptive rule raises
  dps, as it did up to 331 on P4 k = 40.
- **The P4 margin at k = 2** (1 %) is decided exactly, so no float can flip it.

## Time
- Build and smoke: ~1.5 h.
- Run: ~30 min wall.
- Aggregate, checks and thresholds: ~20 min.
- Macros, compilation and the documentation fix: ~30 min.
- Stage 4 (README, prose-author, gate): ~1 h.
- Closing (gates, history, panel): ~20 min.
