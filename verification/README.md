# verification/

Independent re-derivations, written from the definitions and displayed equations of the
manuscript. They import nothing from `experiments/` and use numpy only.

| Script | Checks |
|---|---|
| `build_from_scratch.py` | shared primitives: graphs, Hamiltonians, distances, spectral projectors, separation parameters (not a test) |
| `verify_examples.py` | Examples 7.1–7.3: every displayed parameter and every table row, against the true `S`, `F̄` and a sampled `F_max` in double precision. Called by `experiments/2026-09-18_interval_certification/code/crosscheck.py`. Ends with `ALL EXAMPLE CLAIMS REPRODUCED` |
| `verify_prop_comparacao_round3.py` | Section 8.2, Proposition 8.1 (comparison with a bound built from a global spectral enclosure, on the P4 family of Example 7.3): the constants `a⋆` and `b_k`, the lower bound of the proposition for every admissible triple, and the printed values 1.444…, 0.945…, 1.074…, 0.410… and 2.585·10⁻⁹. Ends with `NO VIOLATION` |

```sh
cd verification
python verify_examples.py
python verify_prop_comparacao_round3.py
```

The docstrings follow the authors' working drafts, not the printed manuscript:
"Section 9" and "Example 9.x" in `verify_examples.py` are Section 7 and Examples 7.x; the
equation labels (8.2), (8.3) and the line numbers in `verify_prop_comparacao_round3.py` refer to
the working source. The remark in `build_from_scratch.py` that the paper ships no code predates
`experiments/`. Double-precision values far below machine precision relative to `‖H‖` (for
example the true `S` of the P4 k = 10 row, printed as 2.1·10⁻³⁴) are round-off; the resolved
multiprecision values are in `experiments/2026-09-18_exact_benchmarks/results/aggregate.csv`.
