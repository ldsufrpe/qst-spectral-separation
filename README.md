# qst-spectral-separation

Code, aggregated results and run logs for the computations in the manuscript
*"Localization of Spectral Projectors from Block Data, with an Application to Quantum State
Transfer"* by Leon D. da Silva, Eberson Ferreira da Silva and Tiago Mendonça Lucena de Veras
(Departamento de Matemática, Universidade Federal Rural de Pernambuco). The manuscript is in
preparation and is not part of this repository.

## What this is

The paper considers a real symmetric matrix `H` supported on a graph, split as
`H = H_in + H_out` with `H_in` block diagonal along a partition of the vertices,
`‖H_out‖ ≤ δ`, and block spectra separated by `g > 2δ`. It proves that the column norms of the
spectral projector `P_i` attached to each block decay exponentially in the graph distance to
that block, with constants computed from block data alone: no eigendecomposition of `H` and no
enclosure of its global spectrum. Summing over blocks bounds the off-diagonal entries of
`e^{-iτH}` uniformly in `τ ≥ 0`, and hence the fidelity of single-excitation quantum state
transfer.

The theorems are analytic. This repository holds the computations that go with them:

1. a rigorous certification of every number printed in the three example tables of Section 7
   (Examples 7.1–7.3);
2. supporting computations that compare the certificates with exact values on the same
   families, and apply them to two further families (a binary tree of modules and a
   mass–spring stiffness matrix);
3. two independent re-derivations in plain double precision (`verification/`).

It is a curated export of the authors' working directory, limited to what is needed to rerun
these computations.

## Findings

- **Tables of Section 7.** All 42 printed entries (14 rows × `Σ1`, `F_esp`, `B²`) are valid
  upper bounds and the tightest upward roundings to four significant figures. The 17
  parameters displayed with "≈" round correctly, and each tabulated `κ` lies in the argmin of
  `Σ1` on the decimal grid of step 10⁻⁴. Arb ball arithmetic (python-flint 0.9.0, FLINT 3.6.0)
  at 256 bits, with identical verdicts at 512 bits.
- **Certificates against exact values.** On all 14 instances every certified inequality holds
  against multiprecision diagonalization, block masses included. The certificates are
  conservative: `Σ1/S` ranges from 4.3·10⁵ to 3.0·10²⁵⁹, because the true decay rate is 4–5
  times the certified rate `γ_κ` on the K4 and dimer families and about 19 times on the P4
  family.
- **Binary tree of K4 modules.** The radial growth condition (GR) holds with `a = ln2/2` and
  `C = 3√2` uniformly in depth, so the prefactor of Theorem C does not depend on the number of
  modules. For a fixed pair, the diameter bound overtakes the Theorem C bound between depths 13
  and 14.
- **Mass–spring stiffness matrix.** The certificates bound the time-uniform displacement
  transmissibility `sup_τ |cos(τ√K)_{ts}|` without modal analysis; they hold for every size
  tested, with certificate-to-exact ratios from 10⁷ to 10⁴⁷.

Each experiment's `README.md` gives the full account, including the hypotheses that were
refuted.

## Repository structure

```
.
├── README.md, LICENSE, LICENSE-DATA, CITATION.cff, requirements.txt
├── reproduction_report.md          clean-environment rerun of every command below
├── experiments/
│   ├── 2026-09-18_interval_certification/   tables of Section 7, Arb certification
│   │   └── data/examples_sec7.tex           Examples 7.1–7.3 (with tables), verbatim extract
│   ├── 2026-09-18_exact_benchmarks/         certificates vs exact values; draft figure
│   ├── 2026-09-18_binary_tree_family/       (GR) constants and Theorem C on a tree
│   └── 2026-09-18_mode_localization/        mass–spring stiffness matrix
│       (each: README.md, plan.md, experiment_log.json, env/, code/, results/)
└── verification/                    independent double-precision re-derivations
```

## Requirements

Python 3.12 (results produced with 3.12.3 on Linux). A LaTeX installation is needed only to
regenerate the figure (the SciencePlots `science` style sets `text.usetex`).

```sh
python3 -m venv .venv
. .venv/bin/activate
pip install -r requirements.txt
```

## Reproducing the results

Run `interval_certification` first: `exact_benchmarks` and `binary_tree_family` read its
`results/aggregate.csv`. Commands are run from `experiments/<id>/code/` unless stated.

| Paper item | Experiment | Commands | Produces | Time |
|---|---|---|---|---|
| Tables of Section 7, displayed parameters, choice of `κ` | `interval_certification` | `python run.py --prec 256`, `python run.py --prec 512`, `python crosscheck.py` | `results/aggregate.csv`, `params.csv`, `kappa_audit.csv` (and `_prec512`) | 3 min |
| Certificates vs exact values on Examples 7.1–7.3 | `exact_benchmarks` | `python run_benchmark.py`, `python analysis.py` | `results/aggregate.csv`, `masses.csv`, `decay_fits.csv` | 14 min |
| Figure `itf_vs_bounds` (draft, not yet in the manuscript) | `exact_benchmarks` | `python run_sweep.py`, `python make_figure.py` | `results/sweep.csv`, `figures/itf_vs_bounds.pdf` | 3 min |
| Binary tree of K4 modules | `binary_tree_family` | `python run.py --rmax 12 --rexact 5`, `python fixed_pair.py` | `results/aggregate.csv`, `fixed_pair.csv`, `radial_constants.csv`, `theoremC_choice.csv` | 28 min |
| Mass–spring stiffness matrix | `mode_localization` | `python run.py --kmax 10` | `results/aggregate.csv`, `participation.csv` | 1 min |
| Proposition 8.1 (Section 8.2) | — | from `verification/`: `python verify_prop_comparacao_round3.py` | stdout, ends with `NO VIOLATION` | 14 s |

All of these commands were rerun from a clean export of this tree in a fresh virtual
environment; `reproduction_report.md` has the exit codes, the wall times above (AMD Ryzen 7
3700U laptop, with the last four experiments running in parallel) and the file-by-file
comparison with the versioned results.

## About the records

- `plan.md` is the plan of each experiment, written before it ran; its hypotheses are checked
  in the experiment's `README.md`, and some were refuted.
- `experiment_log.json` is append-only. Each run records its parameters, input and output
  hashes, wall time, peak memory and status. The fields `provenance_note` and `code_commit`
  refer to the authors' internal tracking (review-item identifiers such as R-C4, commits of the
  private working repository) and are kept verbatim.
- Running a script **modifies versioned files**: it appends a run to `experiment_log.json` and
  rewrites `env/requirements.txt` and `env/system_info.json` with your environment. Work on a
  copy, or use `git restore` afterwards, to keep the recorded state.
- CSV files are stored byte-exact (`.gitattributes`), so the SHA-256 values in the logs can be
  checked against them.
- All computations are deterministic; there are no seeds.
- `families.py`, `bounds.py`, `spectral_mp.py` and `logutil.py` are duplicated across the
  experiments on purpose (each run logs the hash of its own copy); the copies are identical.
- `interval_certification/data/examples_sec7.tex` is a verbatim extract of Examples 7.1–7.3 of
  the manuscript source, tables included; `run.py` parses the printed values from it.

## Citation

See `CITATION.cff`. The manuscript is in preparation; cite this repository together with the
paper once it appears.

## License

Code: MIT (`LICENSE`). Data, results and figures (CSV files, experiment logs, figures):
Creative Commons Attribution 4.0 International (`LICENSE-DATA`).
