# Examples

Runnable scripts, from the repository root (`python examples/<script>.py`). Outputs go to
`results/examples/`.

| Script | Shows |
|---|---|
| [`01_baseline_projection.py`](01_baseline_projection.py) | no-screening projection of India; every compartment of the national model |
| [`02_screening_coverage.py`](02_screening_coverage.py) | HPV-DNA testing at 0–100 % coverage in India, with 95 % intervals, the elimination year and the paper's Fig 3 |
| [`03_compare_screening_tests.py`](03_compare_screening_tests.py) | Pap, VIA, LBC and HPV-DNA testing at equal coverage; the paper's Fig S3 rebuilt; a custom sensitivity |
| [`04_change_parameters.py`](04_change_parameters.py) | what-if runs: changed transmission, progression, start year, a parameter file of your own; effect on the fit |
| [`05_low_level_api.py`](05_low_level_api.py) | model, scenarios, draw schemes and ensembles without `project()` |
| [`06_compare_states.py`](06_compare_states.py) | all 36 state/UT models side by side |
| [`07_model_fit.py`](07_model_fit.py) | fitted models vs the observed registry rates, 2012–2019 (Fig 1, Fig 2) |

Every figure is drawn in the paper's style (same fonts, sizes, colours and legends).
