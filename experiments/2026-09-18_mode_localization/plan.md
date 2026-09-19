# Plan — 2026-09-18_mode_localization

**Date:** 2026-09-18 · **Scale:** canonical · **Mode:** numerical
**Parent:** opportunities.md #1 (LAA framing) · **Delegated finding:** R-B8
**Supersedes:** `plan_mini.md` (approved 2026-09-18)

## Question
For the stiffness matrix of weakly coupled substructures on graded elastic supports,
do the block-projector bounds certify two quantities a priori, with no modal analysis?
(i) The modal participation `‖P_i e_s‖²` of each degree of freedom in each
substructure's band of ω². (ii) The time-uniform displacement transmissibility
`sup_τ |[cos(τ√K)]_{ts}|` between distant substructures. How close do the
certificates come to the exact modal values?

## The linear-algebra fact that carries the example
For real symmetric A with distinct eigenvalues λ and spectral projectors Π_λ,
`S(s,t) = Σ_λ |(Π_λ)_{ts}| = max_{‖f‖∞ ≤ 1} |f(A)_{ts}|`.
The bound ≤ is the triangle inequality, word for word as in Prop. `reducao`. Equality
comes from choosing `f(λ) = sign((Π_λ)_{ts})`. Theorem B therefore bounds the (t,s)
entry of **every** bounded function of A, not just `e^{−iτA}`: `cos(τ√K)` for
displacement response, `e^{−iτK}`, spectral filters, and so on. This is a one-line
remark for §2, and it goes to the next address-review wave (it is not an edit here).
The experiment checks it numerically: the sampled cos-transmissibility never exceeds
`S`, and the sign-function f attains `S` to 1e-50.

## Construction (fixed here, can be revisited at the Stage 2 smoke test)
- Unit masses, M = I. For a diagonal M the pattern survives in `M^{−1/2}KM^{−1/2}`,
  which the README notes.
- There are k substructures, each a chain of 4 masses with internal springs
  `k_in = 1`. Substructure i sits on grounding springs `κ_i = 1 + 6i` (graded support),
  and consecutive substructures are joined end to end by a coupling spring `c = 0.1`.
- `K = Σ springs`, with the diagonal equal to the sum of the attached stiffnesses (the
  coupling spring adds c to both end diagonals). K is supported on a path of 4k
  vertices, with d = 2 and μ = max(k_in, c) = 1.
- `K_i` = path Laplacian of P4 + `κ_i I` + c at the bridge ends. Its width is about
  3.41 + c and its closed-form spectrum is checked, which gives g ≈ 6 − 3.41 − c ≈ 2.49
  and δ = c = 0.1. So (H2) holds with a large margin. g is not tuned to be marginal.
- s is the free end of substructure 0 and t is the free end of substructure k−1, for
  k = 3…10.
- **Physical reading:** a structure on a graded elastic foundation with weak
  inter-segment joints. Each ω²-band belongs to one segment, and the certificate
  bounds how much of an initial displacement at s is ever seen at t, and how much each
  DOF participates in a foreign band, without computing a single mode.

## Hypothesis
Supported: the hypotheses hold, `‖P_i e_s‖ ≤ L_i(s)` everywhere,
`sup_τ |cos(τ√K)_{ts}| ≤ S ≤ Σ1`, and `Σ1/S` is of the same order as in the
quantum-chain families.
Refuted, or fallback: (H2) needs contrived stiffness ratios. If a support step of 6
(relative to internal springs of 1) is judged unrealistic at review, the example is
reported as a structured-matrix illustration without the physical claims, or R-B8
falls back to the f(A)-decay framing.

## Methodological critique of the plan-mini
1. **Is this just Example 7.3 renamed?** Both have P4 modules. The diagonal here
   differs (the Laplacian structure gives endpoint diagonals 1 + κ_i and interior ones
   2 + κ_i, not constant), and so does the observable (cos(τ√K) instead of e^{−iτH}).
   To avoid "same matrix, new story", the example also reports the band participation
   `‖P_i e_s‖²`, a quantity with its own engineering meaning (the modal participation
   factor of the band), which the quantum examples never print.
2. **Is ω² the right variable?** The hypotheses apply to the spectrum of K, i.e. ω²,
   and the bands in ω follow by monotonicity. Stated explicitly.

## Code reuse
`spectral_mp.py` copied from `exact_benchmarks` (SHA logged). New: `structure.py`
builds K, and `run.py`.

## Verification order
1. K is SPD, the pattern matches the path graph, d = 2 and μ = 1.
2. The exact block spectra and (H1), (H2).
3. Lemma aglo counts.
4. `max_f` equality via the sign function.
5. `‖P_i e_x‖ ≤ L_i(x)` and `cos-sup ≤ S ≤ Σ1`, with a stop on violation.
6. The ratios.

## Files
- **To create:** `code/*`, `results/{aggregate.csv, participation.csv, viz_schema.json}`,
  `experiment_log.json`, `README.md`, `env/*`.

## Time
Build ≈ 1 h, run < 5 min (n ≤ 40 in mpmath).
