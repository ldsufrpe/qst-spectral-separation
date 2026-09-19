# Reportable narrative

On the fourteen instances of Examples 7.1–7.3, the certificates Σ₁, F_esp and B² bound the exact quantities they control but sit far more than a few orders of magnitude above them, so the hypothesis of near-tightness is refuted. The exact ITF S², time average and block-projector masses, together with a sampled supremum of the transfer probability, come from multiprecision diagonalization at 60 to 331 digits with checked convergence and serve only as calibration, since no certificate uses them. No certified inequality fails on any instance, block masses included. The sampled supremum attains the ITF to three significant figures on the K4 and dimer instances and reaches 0.97 of it on P4, so the gap lies between certificate and ITF, where Σ₁/S grows from about 4·10⁵ (K4, k=3) to 3·10²⁵⁹ (P4, k=40). The fitted decay rate of transport, 5.4 to 5.9 per unit distance for K4 and the dimer and 4.4 for P4, exceeds the certified γ_κ of 0.8 to 1.4 and 0.22 to 0.24 by a factor of four to five and of about nineteen. The certified rate is normalized by the internal coupling and sees the weak bridges only through η, the floor at which L_i stays on nearby blocks, so the certified block-mass profile has the true shape without the true rate. A Stark-ladder factorial correction halves the fit residual, which with 3 to 10 points per family is suggestive only.

# Operational log

**Experiment:** `2026-09-18_exact_benchmarks` · mode numerical
**Pre-registered plan:** `plan.md` (2026-09-18; its "Tables 3–5" are the three tables of
Section 7 under an earlier numbering)

## What was computed
For each of the 14 rows of Table 2 of Section 7 (Examples 7.1–7.3, printed as three tables when this ran; s and t as in §7,
κ = κ_tab), the following exact
quantities were computed by multiprecision diagonalization (mpmath `eigsy`):
- `S(s,t)` and ITF `S²`, grouped by distinct eigenvalue with tolerance 10⁻⁴⁰
- `F̄ = Σ c_λ²`
- a sampled lower estimate `F_samp` of `F_max`, with `F_max ∈ [F_samp, S²]`
- the block masses `‖P_i e_s‖` and `‖P_i e_t‖` for every block (`results/masses.csv`)

The certificates are the Arb values of `../2026-09-18_interval_certification`, and the
mpmath re-evaluation agrees with them to ≤ 3.6·10⁻¹⁶ relative. The diagonalization is
calibration only. No certificate uses it.

**Precision:** the working precision is adaptive. If `S < 10^−(dps−30)`, the run is repeated
with more digits, and a final run at +20 digits must agree on S to within 10⁻⁶ relative. The
dps actually used ranged from 60 to 331; the P4 k=40 row has `S = 4.3·10⁻²⁷¹`. Where
`S > 10⁻¹⁰`, float64 `eigh` agrees to ≤ 1.6·10⁻¹⁰. The float64 value of S for P4 k=10
printed by `../../verification/verify_examples.py` (2.1·10⁻³⁴) is round-off noise; the
resolved value is 9.8·10⁻⁴⁴. That changes no verdict.

## Outcome
| family | k | D | exact ITF `S²` | `F_samp` / ITF | printed `B²` | `Σ1/S` | `B²/S²` |
|---|---|---|---|---|---|---|---|
| K4 | 3 | 5 | 5.48e-16 | 1.00 | 1.01e-4 | 4.3e5 | 1.8e11 |
| K4 | 5 | 9 | 3.30e-33 | 1.00 | 1.29e-5 | 6.2e13 | 3.9e27 |
| K4 | 10 | 19 | 7.63e-82 | 1.00 | 7.75e-13 | 3.2e34 | 1.0e69 |
| dimer | 3 | 5 | 8.03e-18 | 1.00 | 5.29e-5 | 2.6e6 | 6.6e12 |
| dimer | 12 | 23 | 6.12e-110 | 1.00 | 8.45e-21 | 3.7e44 | 1.4e89 |
| P4 | 10 | 39 | 9.63e-87 | 0.97 | 6.12e-2 | 2.5e42 | 6.4e84 |
| P4 | 40 | 159 | 1.87e-541 | 0.97 | 1.66e-22 | 3.0e259 | 8.9e518 |

All 14 rows are in `results/aggregate.csv`. Every structural and inequality check passes on
every row:
- Lemma `aglo` eigenvalue counts per `I_i`
- `S ≤ Σ1`, `S² ≤ B²`, `F̄ ≤ Σ2` and `F_samp ≤ S²`
- `‖P_i e_x‖ ≤ L_i(x)` for every block
- agreement of Σ1 with the Arb values
- precision convergence

The sampled supremum attains the ITF to 3 s.f. on the K4 and dimer rows (0.97 on P4). The ITF
is therefore essentially the attained maximum fidelity here, and the gap sits between the
certificate and the ITF, not between the ITF and `F_max`.

**Block masses.** For K4 k=5, `‖P_i e_s‖` falls from 2.0·10⁻⁴ (block 1) to 2.1·10⁻¹⁸
(block 4), while `L_i(s)` stays at its floor `η = 5.0·10⁻³` on blocks 1–3 and only then
decays (1.8·10⁻³ on block 4). The mass ratio `L_i/‖P_i e_s‖` grows from 26 to 8.5·10¹⁴. The
profile has the right shape (monotone, maximal on the home block) but not the right rate.

## Decay rates (`results/decay_fits.csv`, `code/analysis.py`)
| series | k range | fitted `γ_true` per unit D | certified `γ_κ` | RMS residual, exponential | RMS residual, Stark ladder |
|---|---|---|---|---|---|
| K4 (table rows) | 3–10 | 5.43 | 0.81–1.15 | 0.83 | 0.39 |
| K4 (sweep a) | 3–12 | 5.57 | 0.81–1.16 | 1.08 | 0.51 |
| dimer | 3–12 | 5.91 | 1.00–1.41 | 1.31 | 0.63 |
| P4 | 10–40 | 4.40 | 0.22–0.24 | 8.61 | 6.42 |

The ratio `Σ1/S` grows because the true decay rate is 4–5× the certified `γ_κ` in the K4 and
dimer families and about 19× in P4. The certified rate `γ_κ = log(1 + κĝ)` is set by
`ĝ = (g − 2δ)/(2dμ)`, which is normalized by the *internal* coupling `dμ`. The bridge weight
enters only through `η`, which is the floor of `L_i`. A two-parameter Stark-ladder model
(`ln S = c0 + (k−1) ln q − ln (k−1)!`, the leading perturbative form when the site energies
grow linearly along the chain) halves the residual in K4 and dimer. The data therefore show
faster-than-exponential curvature on top of the rate gap. With 3–10 points per family, this
is suggestive, not decisive.

## Figure 1 of the manuscript: `figures/itf_vs_bounds.pdf`
Panels (a), (b) come from the sweep in `results/sweep.csv`. Panels (c), (d) plot the binary tree
of modules from `../2026-09-18_binary_tree_family/results/aggregate.csv` (root to a leaf at depth
r = 1…12) and `fixed_pair.csv` (a fixed pair at D = 17 while the tree grows to k = 131 071
modules); no computation is repeated. K4 family, w = 0.1. (a) Δ_mod = 24, k = 3…12: the ITF falls from 5.5·10⁻¹⁶ to 1.0·10⁻¹⁰²,
while `min{1,Σ1²}` falls from 1.0·10⁻⁴ to 2.3·10⁻¹⁶. (b) k = 5, `g/(2δ)` from 100 down to
1.05: `Σ1²` and `b_η²` reach 1 at `g/(2δ) ≤ 2` and ≤ 1.2 respectively, and `B_E² ≈ 0.54` is
the only informative bound near the edge of (H2). Meanwhile the ITF stays at 8·10⁻¹⁹. κ is
grid-optimal (step 10⁻³, 1/2 on plateaus). Element list: `figures/figure_manifest.json`.
The rendered PDF/PNG is not versioned in this release; `code/make_figure.py` regenerates it
from those CSV files and needs a LaTeX installation (the SciencePlots `science` style sets `text.usetex`).

## Hypothesis check (plan.md)
**Refuted as stated:** "Σ1/S stays within a few orders of magnitude". The ratio ranges from
4.3·10⁵ (K4 k=3, the smallest instance) to 3.0·10²⁵⁹ (P4 k=40). It grows exponentially in D
at a rate close to `γ_true − γ_κ`. The certificates are valid everywhere, and the spatial
profile of `L_i` tracks the shape of the true block masses, but the certified rate is 4–19×
too slow on these families.

## Reproduce
```bash
cd experiments/2026-09-18_exact_benchmarks/code
python run_benchmark.py      # ≈ 15 min (P4 k=40 ≈ 10 min); writes aggregate.csv, masses.csv
python run_sweep.py          # ≈ 3 min
python analysis.py; python make_figure.py
```
Run `../2026-09-18_interval_certification` first: `run_benchmark.py` reads its
`results/aggregate.csv`.
Deterministic, so there are no seeds. Environment: `env/requirements.txt`,
`env/system_info.json`.

## Log notes (`experiment_log.json`, 51 runs, 0 failures)
- r0001 is the smoke run.
- r0002–r0006 and r0010–r0018 are the first full pass. Its CSV writer crashed on
  heterogeneous columns after all runs had been logged. The writer was fixed and the pass
  repeated as `--tag final`, which produced the current `aggregate.csv` and `masses.csv`.
- r0007–r0009 are three points of a first sweep launch. It was stopped to add file locking to
  the log helper before two processes wrote to the same log, and it was superseded by the full
  sweep.

## Pointers
- Code: `code/run_benchmark.py`, `code/spectral_mp.py`, `code/bounds.py`, `code/run_sweep.py`,
  `code/analysis.py`, `code/make_figure.py`, `code/families.py`, `code/logutil.py`.
- Results: `results/aggregate.csv`, `results/masses.csv`, `results/sweep.csv`,
  `results/decay_fits.csv`. Raw stdout logs are not versioned.
