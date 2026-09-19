# Reportable narrative

On a binary tree of K4 modules, where the community count grows exponentially with distance and the chain corollary does not apply, (GR) holds with a = ln2/2 and C = 3√2 uniformly in depth, so the Theorem C prefactor is independent of the number k of modules, while the diameter-corollary prefactor C1^unif grows linearly in k. The constants come from counting modules per tree distance, checked by exact all-sources breadth-first counts up to depth 12 (32 764 vertices) and by the exact multiprecision ITF to depth 5. The counted C(ln2/2) rises from 1.5 to 2.969, never above 3√2, the exact C_θ reaches at most 0.42 of its (GR) bound, and with κ = 0.85 and β = 0.73 the fixed prefactor 2.1·10⁸ makes the Theorem C bound informative from D = 15 (6.5·10⁻²) down to 3.0·10⁻⁸ at D = 25. The exact ITF lies below every certificate. Between root and deepest leaf the diameter bound undercuts the Theorem C bound at every depth, because the module count grows only like e^{(ln2/2)D} and the full rate absorbs it, so there Theorem C gains k-independence and no smaller number. For a fixed pair at D = 17 in trees up to 131 071 modules, the diameter bound grows like k² and crosses the Theorem C value 3.49·10⁻³ between depths 13 and 14, the regime where uniformity pays. The sharp Σ₁ is the smallest certificate throughout and, as in the chain families, sits orders of magnitude above the exact ITF.

# Operational log

**Experiment:** `2026-09-18_binary_tree_family` · mode numerical + exact combinatorics
**Pre-registered plan:** `plan.md` (2026-09-18). One correction at the smoke test: `D = 2r+1`, not `2r+2`.

## Family
K4 modules (coupling 1) sit at the nodes of a complete binary tree of depth r, so there are
k = 2^{r+1} − 1 modules and n = 4k vertices. Module ℓ (BFS label) has site energy 24ℓ.
Bridges of weight w = 0.1 run from parent vertex 1 or 2 to child vertex 0; they form a
matching. The hypotheses and parameters, all checked in code:
- (H1): `δ = ‖H_out‖ = 0.1`, since the bridges form a matching.
- (H2): `g = 20`, `W = 4`, from the same blocks as Example 7.1.
- (H3): `D0 = 1`, `M0 = 4`.
- `d = 4`, `μ = 1`, `η = 5.03·10⁻³`, `ĝ = 2.475`, `ϱ = 1.1350`.

s is vertex 3 of the root and t is vertex 3 of the leftmost leaf, with `D = 2r + 1`
(checked by BFS). No chain indexing of the modules exists, so Corollary `cadeiasradiais`
does not apply.

## (GR) constants
- **Analytic, used in Theorem C.** Two facts hold for x in module j:
  - `dist(x, V_i) ≥ 2 d_T(j,i) − 1`, because each module transition after the first costs a
    bridge plus an internal hop.
  - A binary tree has at most `3·2^{ρ−1}` modules at tree distance ρ.

  Hence `N_x(R) ≤ 3·2^{(R+1)/2} − 2 ≤ 3√2·2^{R/2}`, so (GR) holds with `a = ln2/2 ≈ 0.3466`
  and `C = 3√2 ≈ 4.243`, uniformly in r. This bound is not stated in the current
  manuscript; the two facts above are its proof.
- **Exact, by all-sources BFS, r = 1…12 (n up to 32 764).** `C(ln2/2) = max_{x,R} N_x(R)e^{−aR}`
  takes the values 1.5, 2.0, 2.25, 2.5, 2.625, 2.75, 2.8125, 2.875, 2.906, 2.938, 2.953,
  2.969 (increasing, apparently toward 3). It never exceeds 3√2, which the code asserts for
  every r and every a ≥ ln2/2. The sup over R is taken over the finite set of jump points, so
  it is exact, not sampled.
- **Check of the theorem's inequality** `C_θ ≤ Cθ/(θ−a)`: over 1512 pairs (r, θ > a) the
  exact `C_θ` reaches at most 0.42 of the bound (r = 12, θ = 0.48).

## Theorem C constant
The parameters minimize the distance at which `𝒞e^{−2βD} < 10⁻²`, over κ ∈ {0.50,…,0.95},
a ≥ ln2/2 and β on a 0.005 grid:
- `κ = 0.85`, `β = 0.73`, `a = ln2/2`, `α = γ_κ − β = 0.4026 > a`
- `γ_κ = 1.1326`, `A0 = 15.13`
- `𝒞 = A0⁴e^{2βD0}(Cα/(α−a))² = 2.10·10⁸`, and `D* = 16.3`

With the exact family `C_α` (10.93, the sup over r ≤ 12) in place of the (GR) bound, the
same inequality (eq:radialF) gives `2.70·10⁷`. That value is reported as observed, not as
proven uniform.

## Certificates vs exact (`results/aggregate.csv`)
| r | k | D | `Σ1²` (sharp) | Cor. diametros with `C1^unif` | Theorem C, `𝒞e^{−2βD}` | exact ITF |
|---|---|---|---|---|---|---|
| 1 | 3 | 3 | 1.02e-4 | 1 | 1 | 3.19e-8 |
| 3 | 15 | 7 | 1.04e-4 | 1 | 1 | 8.66e-28 |
| 5 | 63 | 11 | 2.02e-7 | 3.92e-2 | 1 | 2.85e-53 |
| 7 | 255 | 15 | 5.95e-10 | 6.16e-5 | 6.47e-2 | — |
| 9 | 1023 | 19 | 8.03e-13 | 1.08e-7 | 1.88e-4 | — |
| 12 | 8191 | 25 | 3.24e-18 | 6.63e-12 | 2.95e-8 | — |

`C1^unif` (eq:C1unif) grows with k, from 20 (r = 1) to 1.26·10⁷ (r = 12), while `𝒞` is
fixed. All checks pass wherever the exact ITF exists (r ≤ 5, dps 60, run by the shared
benchmark code): `ITF ≤ Σ1² ≤ C1^unif bound` and `ITF ≤ Theorem C bound`. Beyond r = 5, n
exceeds the multiprecision budget (≈ 50 min at r = 6), so only the certificates are
reported.

## Hypothesis check (plan.md)
- **Supported:** (GR) holds with `a = ln2/2 > 0` and an explicit C, uniform in r, both
  analytically and in exact counts. An admissible `(κ, β)` with `α0 > a` exists. `𝒞` is
  independent of k, while `C1^unif` grows linearly in k. The Theorem C bound is informative
  from r = 7 (D = 15) and reaches 3·10⁻⁸ at r = 12.
- **Not anticipated.** In this family the non-uniform bound `C1^unif e^{−γD}` is smaller than
  the Theorem C bound at every r. The number of modules grows only like `e^{(ln2/2)D}`, which
  the full rate γ absorbs (`γ − a = 0.79 > β = 0.73`). The advantage of Theorem C here is the
  k-independence of the prefactor, not a smaller number. The sharp `Σ1` is the smallest of the
  three by orders of magnitude. Against the exact ITF, `Σ1²/ITF` ranges from 3.2·10³ (r = 1)
  to 7.1·10⁴⁵ (r = 5), the same conservativeness as in `../2026-09-18_exact_benchmarks`.
- **Fixed pair, growing tree** (`code/fixed_pair.py`, `results/fixed_pair.csv`). This regime
  isolates the k-independence. s is the root and t is the leftmost node at depth 8, so D = 17
  stays fixed while the tree grows from r = 8 to r = 16 (k = 511 → 131 071; certificates
  only, BFS from s and t). The Theorem C bound stays at 3.49·10⁻³. The Cor. diametros bound
  grows like k²: 2.5·10⁻⁶, 6.5·10⁻⁴ (r = 12), 1.0·10⁻² (r = 14), 0.17 (r = 16). It crosses
  the Theorem C bound between r = 13 and r = 14. The sharp `Σ1² = 3.19·10⁻¹¹` does not move,
  because the added modules are far from both s and t. This is the regime where Theorem C is
  the better of the two closed-form bounds.
- **Cosmetic:** at r = 1 the κ column reads 0.001 because Σ1 is flat in κ there and the
  float grid search takes the first minimizer. The value of Σ1 is unaffected.

## Reproduce
```bash
cd experiments/2026-09-18_binary_tree_family/code
python run.py --rmax 12 --rexact 5   # ≈ 45 min (r=12 BFS ≈ 15 min, r=5 exact ≈ 9 min)
python fixed_pair.py                 # ≈ 1 min; reads results/theoremC_choice.csv from run.py
```
Run `../2026-09-18_interval_certification` first: `benchlib.py` reads its `results/aggregate.csv`.
Deterministic, so there are no seeds. Environment: `env/requirements.txt`,
`env/system_info.json`.

## Pointers
- Code: `code/tree.py`, `code/radial.py`, `code/run.py`, `code/fixed_pair.py` (≈ 5 min), and copies of `benchlib.py`
  (= `exact_benchmarks/code/run_benchmark.py`), `bounds.py`, `spectral_mp.py`, `families.py`,
  `logutil.py`.
- Results: `results/aggregate.csv`, `results/fixed_pair.csv`, `results/radial_constants.csv` (C(a) and C_θ per r),
  `results/theoremC_choice.csv`. Raw stdout (`results/run.log`) and smoke outputs
  (`--smoke`, prefix `smoke_`) are not versioned.
- Log: `experiment_log.json`.
