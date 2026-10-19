# qst-spectral-separation

Code and data for the numerical examples of the paper

> L. D. da Silva, E. F. da Silva and T. M. L. de Veras,
> *Certified suppression of quantum state transfer between spectrally separated communities*
> (in preparation).

The paper bounds the transfer of a single excitation between communities of a qubit network
whose energy windows are separated. The bounds use only local data of each community and of the
couplings between them. The scripts here evaluate these bounds on concrete networks, certify them
in ball arithmetic (Arb, via python-flint) and compare them with the exact values.

## Contents

| Folder | Content |
|---|---|
| `experiments/2026-10-08_neumann_14rows` | modular chains: Table 2 and Figure 1 |
| `experiments/2026-10-15_rc3_threshold` | chain length from which each bound keeps the fidelity below 2/3 |
| `experiments/2026-10-08_stark_sharpness` | Wannier–Stark ladder: Example 7.1 and Figure 2 |
| `experiments/2026-09-18_exact_benchmarks` | exact reference values used by the experiments above |
| `experiments/2026-09-18_interval_certification` | Combes–Thomas baseline |
| `verification/` | independent check of the baseline in double precision |

Each experiment folder has its own `README.md`, the code (`code/`), the results (`results/`) and
a run log (`experiment_log.json`).

In the code and the CSV files the bounds have short names: `N1` is the bridge bound, `N2` the
walk-sum bound, `N2p_path` the chain bound, `N3x` the contour certificate and `Sigma1` the
Combes–Thomas bound.

## Installation

Python 3.12. LaTeX is needed only to redraw the figures.

```sh
python3 -m venv .venv
. .venv/bin/activate
pip install -r requirements.txt
```

## Running

Run the commands from `experiments/<folder>/code/`. The earlier experiments do not need to be
rerun, because their results are already in the repository.

Table 2 and Figure 1 (about 10 min), in `2026-10-08_neumann_14rows`:

```sh
python run.py --tag final
python aggregate.py --tag final
python analysis.py
python make_figure.py
```

Fidelity threshold (about 15 min), in `2026-10-15_rc3_threshold`:

```sh
python run.py
python aggregate.py --final
python ct_refine.py
```

Example 7.1 and Figure 2 (about 4 min), in `2026-10-08_stark_sharpness`:

```sh
python run.py --part Z
python run.py --part chain
python analysis.py
python make_figure.py
```

In the first two, `run.py` resumes from `results/raw_final.json`; delete it to recompute from
scratch. Running a script appends to
`experiment_log.json` and rewrites the files in `env/`, so use `git restore` afterwards to return
to the recorded state. `reproduction_report.md` reports a clean rerun of these commands.

## License

Code: MIT (`LICENSE`). Data, results and figures: CC BY 4.0 (`LICENSE-DATA`).
To cite, see `CITATION.cff`.
