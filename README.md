# qst-spectral-separation

This repository holds the code, aggregated results and run logs for the computations in the manuscript
*"Certified suppression of quantum state transfer between spectrally separated communities"* by
Leon D. da Silva, Eberson Ferreira da Silva and Tiago Mendonça Lucena de Veras
(Departamento de Matemática, Universidade Federal Rural de Pernambuco). The manuscript is in
preparation and is not part of this repository.

## What this is

The paper studies a single excitation on a network of coupled qubits whose vertices form
communities. The site-energy windows of the communities are separated by `g`, which exceeds
twice the coupling `δ` between communities. The paper then bounds the transfer amplitude
`|⟨t|e^{-iτH}|s⟩|` at every time `τ`. The bounds are computed a priori from the diagonal blocks
of the Hamiltonian and the couplings between them, without its spectrum. The bound on the mass
of each cluster projector decays as `q^N` in the number `N` of bridges, where `q = 2δ/g`. On
chains of communities whose bridges are weaker than a quarter of the separation, it decays
superexponentially. Evaluated in ball arithmetic with the exact
block resolvents, the contour integral behind these bounds gives a computer-assisted
certificate.

The theorems are analytic. This repository holds the computations of the paper's Section 7,
where the bounds are set against exact values:

1. the three families of modular chains (`K4`, dimer and `P4` modules): Table 2, Figure 1 and
   the numbers quoted with them;
2. on the same families, the chain length from which each certificate keeps the fidelity of a
   transferred qubit below the classical value 2/3 (Corollary 6.3), in the paragraph after
   Table 2 and the comparison of decay rates;
3. the Wannier–Stark ladder (Example 7.1): Figure 2 and the numbers quoted with it.

The repository is a curated export of the authors' working directory. It is limited to what
these computations need.

## Findings

- **Modular chains (14 instances, Table 2, Figure 1).** Every certificate is an Arb upper bound
  computed from block data alone. The exact overlap sum `S(s,t)` comes from multiprecision
  diagonalization and is used only for calibration. The findings:
  - no certified inequality fails on any instance, for the total sum or for any block mass;
  - the contour certificate lies 0.6 to 5.5 orders of magnitude above `S` and captures 0.98 to
    0.99 of the exact decay rate in the number of modules;
  - the Combes–Thomas baseline lies 5.6 to 259.5 orders above `S` and captures 0.05 to 0.17 of
    that rate.

  The rate fractions are secants between the shortest and the longest chain of each family.
- **Fidelity threshold (121 chains, `k = 2` to 10).** By Corollary 6.3, a certificate below
  `√2 − 1` keeps the fidelity of a transferred qubit below 2/3 at every time. Each verdict
  `(Σ₁ + 1)² < 2` is decided in exact rational arithmetic on the certified upper endpoints.
  - On the `P4` chains of Table 1, the contour certificate, the walk-sum and chain bounds and
    `S` lie below the threshold at every `k`. The bridge bound exceeds it at `k = 3`, and the
    Combes–Thomas bound at `k = 3` to 8.
  - On the `K4` and dimer chains of Table 1, `η ≤ 0.005` caps every certificate. The largest
    `Σ₁` on the scanned chains (`k = 2` and 3) is 0.010, and the longer chains of Table 2 also
    lie below the threshold.
  - With the bridges strengthened to `q = 0.40` (`w = 4` on `K4`, `w = 2.76` on the dimer),
    only the contour certificate stays below the threshold at `k = 2`, as `S` does. The
    Combes–Thomas bound, with `μ = w` in its rate, exceeds it at every `k ≤ 6` on the dimer
    chains and every `k ≤ 10` on the `K4` chains.

  Every certificate lies above `S` on all 121 chains. On the three chains shared with Table 2,
  the values agree with those of the other experiments in 18 of 18 comparisons
  (`results/reproduction.csv`). The files also cover the rest of the grid `q = 0.20` to 0.45,
  which the paper does not quote.
- **Wannier–Stark ladder (Example 7.1, Figure 2).** The ladder has site energies `lΔ`, hopping
  `μ` and one site per community; the experiment's files write `F` for `Δ` and `J` for `μ`.
  - On the infinite ladder, the excess of the chain-bound certificate over `S` follows the law
    of Eq. (40) of the paper to within 0.005 in `log10`, for `10 ≤ N ≤ 100` and seven values of
    `Δ/μ` from 5 to 100.
  - On finite chains up to `N = 60`, the walk-sum bound lies 0.9 to 2.4 orders above `S` and
    captures 0.994 to 0.996 of the exact rate.

Each experiment's `README.md` gives the full account, including the hypotheses that were
refuted along the way.

### Names in the files and in the paper

The code and the CSV files keep the internal names of the certificates. The paper calls them:

| In the files | In the paper |
|---|---|
| `N1` | bridge bound (Theorem 5.4) |
| `N2` | walk-sum bound (Lemma 5.5) |
| `N2p_path` | chain bound (Theorem 5.6); `N2p_stated` and `N2p_ladder` are its forms (i) and (ii) |
| `N3x` | contour certificate (Proposition 5.8, Appendix A) |
| `Sigma1` | Combes–Thomas bound (the baseline of Section 7 and Appendix B) |
| `floor` | `Σ_i ‖P_i e_s‖ ‖P_i e_t‖` with the exact masses |
| `N3` | hybrid bound (computed, not in the paper) |

The files of `2026-10-15_rc3_threshold` already use the paper's names in their `certificate`
column: `bridge`, `walk` (walk-sum bound), `chain`, `contour`, `CT` (Combes–Thomas bound) and
`exact` (`S`). In `thresholds.csv`, `k_star` is the smallest scanned `k` from which the
certificate stays below `√2 − 1` at every larger scanned `k`, and `fails` lists the `k` where it
does not. The README and `plan.md` of that experiment cite the paper by its LaTeX labels:

| Label | In the paper |
|---|---|
| `Tab:familias`, `Tab:linhas` | Tables 1 and 2 |
| `lema:aglo` (the `aglo` counts) | Lemma 5.1 |
| `teo:pontes`, `prop:passeios`, `cor:cadeia` | Theorem 5.4, Lemma 5.5, Theorem 5.6 |
| `eq:cadeiamassa` | Eq. (21) |
| `prop:certcontorno` | Proposition 5.8 |
| `teo:principal`, `cor:fent` | Theorem 6.1, Corollary 6.3 |
| `lema:ct` | Lemma B.1 |

## Repository structure

```
.
├── README.md, LICENSE, LICENSE-DATA, CITATION.cff, requirements.txt
├── reproduction_report.md          clean-environment reruns of the commands below
├── experiments/
│   ├── 2026-10-08_neumann_14rows/           Table 2 and Figure 1 (modular chains)
│   ├── 2026-10-15_rc3_threshold/            fidelity threshold 2/3 on the modular chains (Corollary 6.3)
│   ├── 2026-10-08_stark_sharpness/          Example 7.1 and Figure 2 (Wannier–Stark ladder)
│   ├── 2026-09-18_interval_certification/   Combes–Thomas baseline (Sigma1), read by neumann_14rows and rc3_threshold
│   └── 2026-09-18_exact_benchmarks/         exact S and block masses, read by the three later experiments
│       (each: README.md, plan.md, experiment_log.json, env/, code/, results/)
└── verification/                    independent double-precision re-derivation of the baseline
```

## Requirements

Python 3.12 (the results were produced with 3.12.3 on Linux). A LaTeX installation is needed only
to regenerate the figures, because the SciencePlots `science` style sets `text.usetex`.

```sh
python3 -m venv .venv
. .venv/bin/activate
pip install -r requirements.txt
```

## Reproducing the paper's Section 7

The two 2026-10-08 experiments read the committed results of the two 2026-09-18 experiments.
`rc3_threshold` reads those of `neumann_14rows` and of both 2026-09-18 experiments, and imports
their code. None of the experiments it reads needs to be rerun first. Run the commands from
`experiments/<id>/code/`.

| Paper item | Experiment | Commands | Produces | Time |
|---|---|---|---|---|
| Check of the reduced contour solve against the dense one | `neumann_14rows` | `python crosscheck.py --M 300` | `results/crosscheck.csv` | 8 s |
| Table 2, Figure 1, the rate fractions | `neumann_14rows` | `python run.py --tag final`, `python aggregate.py --tag final`, `python analysis.py`, `python make_figure.py` | `results/raw_final.json`, `aggregate.csv`, `blocks.csv`, `rate_shares.csv`, `figures/…_bounds_vs_k.pdf` | 9 min |
| Chain length from which each certificate keeps the fidelity below 2/3 (Section 7, after the rate comparison) | `rc3_threshold` | `python run.py`, `python aggregate.py --final`, `python ct_refine.py` | `results/raw_final.json`, `aggregate.csv`, `thresholds.csv`, `reproduction.csv`, `probe_agreement.csv`, `ct_refine.csv` | 11–18 min |
| Example 7.1, Eq. (40), Figure 2 | `stark_sharpness` | `python run.py --part Z`, `python run.py --part chain`, `python analysis.py`, `python make_figure.py` | `results/z_series.csv`, `chains.csv`, `chain_blocks.csv`, `fits.csv`, `chain_summary.csv`, `rate_shares.csv`, `figures/…_rate.pdf` | 4 min |

The `run.py` of `neumann_14rows` and that of `rc3_threshold` resume from
`results/raw_final.json`. Delete that file first to recompute everything rather than only the
missing tasks. The three `run.py` scripts use 6 worker processes by default (`--workers`).
In `rc3_threshold`, `python run.py --smoke` checks the setup on one chain in a few seconds, and
`aggregate.py --final` also requires all 121 chains and writes its summary into
`experiment_log.json`. `ct_refine.py` re-evaluates the Combes–Thomas failures that lie within
5 % of the threshold on a finer grid of `κ`. It was added after the run, and its README explains
why.

The commands above were rerun from a clean export of this tree in a fresh virtual environment.
`reproduction_report.md` has the exit codes, the wall times (AMD Ryzen 7 3700U laptop) and the
file-by-file comparison with the versioned results.

## About the records

- `plan.md` is the plan of each experiment, written before it ran. The experiment's `README.md`
  checks its hypotheses, and some were refuted.
- `experiment_log.json` is append-only. Each run records its parameters, its input and output
  hashes, the wall time, the peak memory and the status. The fields `provenance_note` and
  `code_commit`, and some file names in `code_files`, refer to the authors' internal tracking:
  - review items;
  - a calibration prototype that is not published;
  - commits of the private working repository.

  They are kept verbatim.
- The published code differs from the code that produced the logged runs in five lines, all
  marked `public export`:
  - both `run.py` files of the 2026-10-08 experiments skip the unpublished calibration
    prototype when they hash their inputs;
  - the `aggregate.py` of `neumann_14rows` and that of `rc3_threshold` read
    `exact_benchmarks/results/aggregate.csv` instead of `final.csv`. The two files are
    byte-identical, and only the former is versioned;
  - `rc3_threshold/exploration/preplan_probe/probe.py` finds the repository root from its own
    location instead of an absolute path on the authors' machine.

  `aggregate.py --oracle` of `neumann_14rows` rebuilds `results/probe_oracle.csv` from an
  internal float probe that is not published. Without the flag, the versioned `probe_oracle.csv`
  is used.
- `rc3_threshold/exploration/preplan_probe/` is the float probe run before that experiment's
  plan, with 300 contour arcs, a `κ` grid of step 1/50 and `S` in double precision. It shaped the
  design and is not evidence. `aggregate.py` reads its output to check that the verdicts agree
  (`results/probe_agreement.csv`). Run from that folder, `python probe.py > probe_output.txt`
  regenerates the output in about 10 minutes.
- The README and `plan.md` of `rc3_threshold` are kept verbatim. They also name files that are
  not published: the superseded first draft of the plan (`plan_mini.md.archived`), the
  manuscript and the script that turns `thresholds.csv` into its macros, the authors' work log,
  and the stdout logs `results/run_*.log`.
- Running a script **modifies versioned files**. It appends a run to `experiment_log.json` and
  rewrites `env/requirements.txt` and `env/system_info.json` with your environment. Work on a
  copy, or run `git restore` afterwards, to keep the recorded state.
- CSV files are stored byte-exact (`.gitattributes`), so the SHA-256 values in the logs can be
  checked against them.
- All computations are deterministic, and there are no seeds.
- `logutil.py` and some numerical helpers are duplicated across the experiments on purpose,
  because each run logs the hash of its own copy.
- The 2026-09-18 experiments were run for an earlier version of the manuscript. Their READMEs,
  docstrings and `data/examples_sec7.tex` use that version's numbering and title. There the 14
  rows are "Table 2" or "Tables 2–4" with Examples 7.2–7.4, and in the 2026-10-08 READMEs they
  are "Table 3". In today's paper they are the rows of Table 2 in Section 7. The figure
  `itf_vs_bounds` of `exact_benchmarks` belongs to that earlier version. Its panels (c) and (d)
  read the binary-tree experiment, which left this tree in the 2026-10 update, together with
  the mass–spring experiment and the check of the earlier Proposition 8.1. All three remain in
  the history at commit `79d4382`.

## Citation

See `CITATION.cff`. The manuscript is in preparation. Once the paper appears, cite this
repository together with it.

## License

The code is under the MIT license (`LICENSE`). The data, results and figures (CSV files,
experiment logs, figures) are under the Creative Commons Attribution 4.0 International license
(`LICENSE-DATA`).
