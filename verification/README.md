# verification/

Independent re-derivations, written from the definitions and displayed equations of the
manuscript. They import nothing from `experiments/` and use numpy only.

| Script | Checks |
|---|---|
| `build_from_scratch.py` | shared primitives: graphs, Hamiltonians, distances, spectral projectors, separation parameters (not a test) |
| `verify_examples.py` | Examples 7.2–7.4: every displayed parameter and every table row, against the true `S`, `F̄` and a sampled `F_max` in double precision. Called by `experiments/2026-09-18_interval_certification/code/crosscheck.py`. Ends with `ALL EXAMPLE CLAIMS REPRODUCED` |

```sh
cd verification
python verify_examples.py
```

The docstrings follow the authors' working drafts, not the printed manuscript:
"Section 9" and "Examples 9.1–9.3" in `verify_examples.py` are Section 7 and Examples 7.2–7.4 of an
earlier version of the manuscript, whose 14 rows are today's Table 2. In today's paper the
certificate checked here (`Σ1`) is the Combes–Thomas baseline of Section 7. The remark in `build_from_scratch.py` that the paper ships no code predates
`experiments/`. Double-precision values far below machine precision relative to `‖H‖` (for
example the true `S` of the P4 k = 10 row, printed as 2.1·10⁻³⁴) are round-off; the resolved
multiprecision values are in `experiments/2026-09-18_exact_benchmarks/results/aggregate.csv`.
