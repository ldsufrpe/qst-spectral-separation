# Reportable narrative

Every number printed in the three example tables of Section 7 is a valid upper bound on the quantity it reports and is its tightest upward rounding to four significant figures, and each tabulated κ minimizes Σ1 on the decimal grid of step 10⁻⁴. Each quantity was recomputed at the printed κ from exact rational inputs in Arb ball arithmetic with outward rounding (python-flint 0.9.0, FLINT 3.6.0) at 256 bits, and every printed value was compared with the exact rational endpoints of its enclosure. All 42 entries of the 14 rows, across the columns Σ1, F_esp and B², are certified upper bounds and tight roundings, with enclosure radii below 10⁻⁷⁵ and identical verdicts at 512 bits. The 17 parameters displayed with ≈ in Examples 7.1–7.3 round correctly, the three printed bounds on b_η² are valid and tight, and every tabulated κ lies in the grid argmin of Σ1. That minimizer is unique in eleven rows. In the other three, Σ1 is constant on the whole grid, so κ = 1/2 is a convention and the bound B² in those rows is set by b_η² instead of Σ1. An independent re-derivation reproduces every claim of the three examples.

# Operational log

**Experiment:** `2026-09-18_interval_certification` · mode symbolic (rigorous arithmetic)
**Pre-registered plan:** `plan.md` (2026-09-18)

## What was certified
- All 14 rows of `Tab:cliques`, `Tab:dimero` and `Tab:p4`, parsed at run time from the
  manuscript text. In this release the parsed file is `data/examples_sec7.tex`, a verbatim
  extract of Examples 7.1–7.3 with their tables. The SHA-256 of the parsed file is in each log
  run. The runs logged on 2026-09-18 parsed the full manuscript source
  (`sha256:c9b80bb6…`), which is the file the extract was cut from (lines 1903–2085). The columns `Σ1`, `F_esp` and `B²`
  were evaluated at the printed `κ_tab`.
- 17 parameters displayed with "≈" in the text of Examples 7.1–7.3, plus 3 printed bounds
  `b_η² ≤ …`.
- The decimal κ search, over 9999 grid points per row (step 10⁻⁴), with Σ1 as the objective.
  §7 does not name the objective; Σ1 is the only printed column that depends on κ.

## Method
- **Library:** python-flint 0.9.0 (FLINT 3.6.0, Arb ball arithmetic). Precision is 256 bits,
  with a 512-bit rerun. Every quantity is a midpoint–radius ball that rigorously contains the
  exact value, and each operation rounds outward.
- **Inputs:** exact rationals (`w`, energies, couplings as `fmpq`). `W_i` and `g` come from
  the exact characteristic polynomial of each block (`fmpz_mat.charpoly().factor()`) and
  FLINT's rigorous root enclosures. No floating-point eigen-solver enters the certified chain.
  Distances are integer BFS values and match `D` in the table.
- **Hypotheses checked in code:** the bridges form a matching, so `‖H_out‖ = w` (H1), and
  `g > 2δ` holds rigorously (H2).
- **Printed-value test:** `ok_upper` means printed ≥ the ball's upper endpoint. `ok_tight`
  means, in addition, that the next smaller 4-significant-figure number already lies below the
  ball's lower endpoint. The comparison uses exact rational endpoints.
- `min` of balls: an ordered or dominated ball is returned exactly. For 18–20 minima over
  overlapping balls, a rigorous hull is used; all of them are equal inter-cluster gaps of
  translated blocks, labelled `gap-min`. The final radii stay below 6.4·10⁻⁷⁶, both absolute
  and relative to 4 s.f.

## Outcome against the hypothesis (plan.md)
| check | outcome |
|---|---|
| printed `Σ1`, `F_esp`, `B²` are upper bounds (14 rows × 3) | 42/42 |
| each printed value is the tightest 4-s.f. upward rounding | 42/42 |
| displayed "≈" parameters round correctly from the ball | 17/17 |
| printed `b_η²` bounds valid and tight | 3/3 |
| `κ_tab` lies in the ball-aware grid argmin of Σ1 | 14/14 |
| 256-bit and 512-bit verdicts identical | yes |
| independent `verification/verify_examples.py` (subprocess) | all example claims reproduced |

The κ audit found a unique grid minimizer in 11 rows. In the other three (K4 k=3,4 and dimer
k=3), Σ1 is flat over the whole grid: every `L_i` sits on its floor `b_i`. `κ = 0.5` is the §7
plateau convention, and in those rows `B²` is set by `b_η²`, not by Σ1 (see `active_bound`).
Hypothesis **supported**. No printed value under-reports a bound.

## Reproduce
```bash
cd experiments/2026-09-18_interval_certification/code
python run.py --prec 256     # ≈ 85 s on the machine below
python run.py --prec 512     # ≈ 100 s
python crosscheck.py         # runs ../../../verification/verify_examples.py
```
Environment: `env/requirements.txt` (pip freeze) and `env/system_info.json` (AMD Ryzen 7
3700U, 8 threads, 9.6 GB, Linux 7.0, Python 3.12.3). The computation is deterministic, so
there are no seeds.

## Files
- `results/aggregate.csv`: per row, the printed value, certified upper endpoint, radius,
  `ceil4`, `ok_upper`, `ok_tight`, active bound and κ audit (256 bits). `aggregate_prec512.csv`
  is the rerun.
- `results/params.csv`: displayed parameters. `results/kappa_audit.csv`: grid argmin per row.
- `results/crosscheck_verify_examples.log`: raw stdout of the cross-check, regenerated by
  `crosscheck.py` and not versioned.
- The smoke run (r0001, `--smoke`: 1 row per family, grid 1/100) writes
  `results/smoke_aggregate.csv`, which is not versioned.
- `data/examples_sec7.tex`: the input parsed by `run.py` (see above).

## Pointers
- Code: `code/families.py`, `code/certify.py`, `code/run.py`, `code/crosscheck.py`,
  `code/logutil.py` (shared helper; identical copies in the sibling experiments).
- Log: `experiment_log.json` (runs r0001 smoke, r0002 256-bit, r0003 512-bit, r0004 crosscheck).
