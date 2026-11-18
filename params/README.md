# Fitted parameters

| File | Rows | Columns |
|---|---|---|
| `states_final.csv` | 36 state/UT models | `State, p, K1, K2, beta_const, beta_Reinfec, gamma, r_I, r_C, mu_C, n_PL, omega, epsilon` |
| `india_final.csv` | national model | the same (without `State`) plus `mu, H_pos` |
| `india_uncertainty_sets.csv` | 500 alternative national parameter sets, best first | `rank, sample_index, beta, beta_reinfection, gamma, r_I, epsilon, mu_C, omega, r_C` |

* Values are stored at full double precision; read them with
  `pandas.read_csv(..., float_precision="round_trip")` (as `hpvcc.parameters` does) to
  reproduce the results exactly.
* State models use μ = 0.01388 and H⁺ = 0.95 (HPV-DNA test).
* The clearance rate r_H and the progression rate θ of screen-detected infection are in
  no file: they are sampled in every simulation, r_H ~ U(r_I, 1) and θ ~ U(0, γ), redrawn
  every year.
* `india_uncertainty_sets.csv` lists the accepted sets of the national Monte Carlo (GLUE)
  analysis; `sample_index` is each set's position among the 10,000 sampled sets and the
  columns not listed keep the values of `india_final.csv`.

Meanings and per-unit values: [`docs/parameters.md`](../docs/parameters.md).

## How the parameters were obtained

Each model was fitted, on the no-screening model, to the crude incidence and mortality
rates of the population-based cancer registry estimates for 2012–2019
([`data/registry/`](../data/registry/)), minimising a weighted relative root-mean-square
error over the two series (details in the paper's Methods). The fitting code is not part of
this repository: every projection runs from the fitted values in this folder.
