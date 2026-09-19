# Reportable narrative

For the stiffness matrix $K$ of weakly coupled mass–spring substructures on a graded elastic foundation, the block-projector certificates bound the band participation of every degree of freedom and the time-uniform displacement transmissibility $\sup_\tau |[\cos(\tau\sqrt{K})]_{ts}|$ without modal analysis, although no fixed factor separates them from the exact values. Chains of four unit masses, with grounding stiffness $1+6i$ on substructure $i$ and coupling springs of stiffness $1/10$, were assembled for $k = 3,\dots,10$ substructures, and hypotheses (H1)–(H2) hold with $g/(2\delta) \approx 13$. The identity $S(s,t) = \max_{\|f\|_\infty \le 1} |f(K)_{ts}|$ turns the certificate into a bound on every bounded function of $K$, the cosine propagator included. Every bound held for every $k$, and the sampled transmissibility reached between 0.82 and 0.94 of $S$, whereas the ratio of the certificate to $S$ grew from about $10^{7}$ at $k=3$ to about $10^{47}$ at $k=10$. The loss is a rate gap: a fit of $\ln S$ against the distance gives a true decay rate of about 3.5 per unit distance, against a certified rate of about 0.41. A factorial correction of the kind expected on the linear ladder of the graded foundation lowers the fit residual only modestly. The experiment therefore supports the validity and the applicability of the certificates to a mechanical operator and does not support their quantitative tightness.

# Operational log

**Experiment:** `2026-09-18_mode_localization` · mode numerical
**Pre-registered plan:** `plan.md` (2026-09-18)

## System
Unit masses, so the eigenproblem is `Kφ = ω²φ`. There are k substructures, each a chain of
4 masses with internal springs `k_in = 1`. Every mass of substructure i sits on a grounding
spring `κ_i = 1 + 6i` (a graded elastic foundation). Coupling springs `c = 1/10` join the
last mass of substructure i to the first of i+1. `K` is supported on the path of `4k`
vertices, with `d = 2`, `μ = 1`, and diagonal entries equal to the attached stiffnesses. The
coupling springs form a matching, so `δ = ‖K_out‖ = c`. The block spectra give
`g = 2.5936` and `max W_i = 3.3989` (the end blocks carry the extra `c` on one diagonal
entry), with `η = 0.0401`, and (H2) holds with `g/(2δ) ≈ 13`. `K` is positive definite, with
`ω²_min = 1.023`. `s` is the free end of substructure 0 and `t` the free end of substructure
k−1, so `D = 4k − 1`.

## The matrix-function statement tested
`S(s,t) = Σ_λ |(Π_λ)_{ts}| = max_{‖f‖∞≤1} |f(K)_{ts}|`. The upper bound is the triangle
inequality (the proof of Prop. `reducao` verbatim). Equality is attained by
`f(λ) = sign((Π_λ)_{ts})`, and the code checks this to working precision. Hence
`sup_τ |[cos(τ√K)]_{ts}|`, the displacement at t after a unit initial displacement at s
released from rest, is bounded by `S ≤ Σ1` for all time.

## Outcome (k = 3…10, multiprecision, dps 60–120 adaptive)
| k | D | exact S | sampled `sup_τ |cos(τ√K)_ts|` | certificate `B` | `Σ1/S` |
|---|---|---|---|---|---|
| 3 | 11 | 1.04e-8 | 9.77e-9 | 8.01e-2 (`b_η`) | 7.9e6 |
| 4 | 15 | 9.63e-14 | 8.86e-14 | 8.01e-2 (`b_η`) | 8.7e11 |
| 5 | 19 | 3.24e-19 | 2.87e-19 | 7.73e-2 | 2.4e17 |
| 6 | 23 | 4.77e-25 | 4.24e-25 | 2.00e-2 | 4.2e22 |
| 7 | 27 | 3.56e-31 | 3.10e-31 | 7.69e-3 | 2.2e28 |
| 8 | 31 | 1.48e-37 | 1.30e-37 | 4.08e-3 | 2.8e34 |
| 9 | 35 | 3.71e-44 | 3.13e-44 | 2.10e-3 | 5.6e40 |
| 10 | 39 | 5.93e-51 | 4.84e-51 | 6.49e-4 | 1.1e47 |

All checks pass in every run: the Lemma `aglo` counts, `S ≤ Σ1`, `F̄ ≤ Σ2`, the masses
`‖P_i e_x‖ ≤ L_i(x)`, the max-f identity, `cos-sup ≤ S`, and precision convergence. The
sampled cos transmissibility reaches 0.82–0.94 of `S`, so `S` is close to the actual
supremum here.

**Band participation** (`results/participation.csv`): `‖P_i e_s‖²`, the fraction of an
impulse at s carried by the modes of band i, against the bound `min{1, L_i(s)²}`. At k=6 the
bound is flat at `η² = 1.61·10⁻³` on bands 1–3, then drops to `4.8·10⁻⁵` on band 5. The
exact values fall from `3.2·10⁻⁸` (band 1) to `1.2·10⁻⁵⁶` (band 5).

## Hypothesis check (plan.md)
- **Supported:** the hypotheses hold with a wide margin (`g/2δ ≈ 13`) and physically ordinary
  parameters. The certificates are valid for this non-quantum operator and bound a quantity
  with a mechanical meaning, with no modal analysis.
- **Not supported as written:** "within a stated factor of the exact modal value". The ratio
  grows from 10⁷ to 10⁴⁷ over k = 3…10, the same pattern as the §7 families
  (`../2026-09-18_exact_benchmarks`). The graded foundation is a linear ladder in the
  diagonal. Prop. `largura` forces the diagonal range to grow at least linearly with the number
  of substructures under (H2), and a monotone ladder is this family's way of meeting it. A fit of
  `ln S` against D gives a true decay rate `γ_true ≈ 3.49` per unit distance, against the
  certified `γ_κ ≈ 0.41` (0.26 on the κ = 1/2 plateau rows k = 3, 4). That rate gap of
  about 8.5× is what drives the ratio. A two-parameter Stark-ladder model,
  `ln S = c0 + (k−1) ln q − ln (k−1)!`, lowers the RMS residual only from 1.52 to 1.11, so
  the data show a mild faster-than-exponential curvature, not a decisive one. The fit uses
  the same code as `../2026-09-18_exact_benchmarks/code/analysis.py`.

## Reproduce
```bash
cd experiments/2026-09-18_mode_localization/code
python run.py --kmax 10      # ≈ 1 min
```
Deterministic, so there are no seeds. Environment: `env/requirements.txt`,
`env/system_info.json`.

## Pointers
- Code: `code/structure.py` (K), `code/run.py`, and copies of `bounds.py`, `spectral_mp.py`,
  `benchlib.py` (= `exact_benchmarks/code/run_benchmark.py`), `families.py`, `logutil.py`.
- Results: `results/aggregate.csv`, `results/participation.csv`. Raw stdout
  (`results/run.log`) and smoke outputs (`--smoke`, `smoke_*.csv`) are not versioned.
- Log: `experiment_log.json` (r0001 smoke, r0002–r0009 full).
