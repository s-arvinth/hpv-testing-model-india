# Code structure

The package `hpvcc` (HPV / cervical cancer) lives in [`src/hpvcc`](../src/hpvcc). The core
modules know nothing about the paper; everything specific to the manuscript lives in the
`hpvcc.paper` subpackage and the top-level [`paper/`](../paper) folder.

```text
src/hpvcc/
├── __init__.py        public API: project, Projection, HPVCervicalCancerModel, ModelParameters, ...
├── constants.py       compartments, years, natural mortality, screening tests, WHO proxy
├── units.py           the 37 units, regions, name maps of the data sources
├── paths.py           repository locations (results/ can be moved with HPVCC_RESULTS)
├── parameters.py      ModelParameters, load_parameters, parameter metadata
├── inputs.py          UnitInputs, RecruitmentInflow (JSON-backed)
├── model.py           ODE right-hand sides + HPVCervicalCancerModel
├── metrics.py         CIR (inflow method), CMR, ensemble interval convention
├── simulation.py      Scenario, DrawScheme family, StochasticEnsemble, EnsembleResult
├── projection.py      project() → Projection (tables, summaries, plots)
├── uncertainty.py     ParameterSetEnsemble → UncertaintyResult
├── observed.py        observed 2012–2019 rates (data/registry), model_fit()
├── data/              population / prevalence sources → data/processed/units/*.json
│   ├── sources.py       readers for data/raw
│   ├── states.py        state/UT inputs
│   └── national.py      national inputs
├── viz/
│   ├── style.py         the paper's figure styles (Fira Sans, SciencePlots base), save, background
│   ├── colors.py        one colour per state, region hue families
│   ├── figures.py       the paper's figure types as functions of data (coverage, fit, tests, uncertainty)
│   └── plots.py         plots of your runs: plot_projection, plot_coverage_bars, plot_fit, plot_compartments
├── paper/             reproduction of the manuscript (hpvcc paper ...)
│   ├── experiments.py   named, cached simulations from paper/experiments.yaml
│   ├── tables.py        master results tables
│   ├── latex_tables.py  the manuscript's LaTeX tables
│   └── figures/         registry of every paper figure; collects the paper's simulations
└── cli.py             the hpvcc command
```

## How a simulation flows

1. **Parameters** — `load_parameters(unit)` reads one row of `params/*.csv` into a frozen
   `ModelParameters` dataclass. Changes produce copies (`params.replace(beta=0.4)`).
2. **Inputs** — `UnitInputs.load(unit)` reads `data/processed/units/<unit>.json`:
   population by age band, rate denominator, the recruitment inflow Λ(t) (a callable) and the
   2002 initial state.
3. **Model** — `HPVCervicalCancerModel(params, inputs, screening_start)` packs the
   parameters into a tuple and integrates `rhs_no_screening` or `rhs_screening` one year
   at a time (`step`, `step_screened`, `simulate`).
4. **Scenarios and draws** — a `Scenario` is a coverage and a sensitivity. A `DrawScheme`
   provides the yearly `(r_H, θ)` of each simulation. `StochasticEnsemble` runs
   scenarios × simulations in a process pool and returns an `EnsembleResult` (arrays of `C`
   and `C_D`, optionally every compartment), which can be saved and reloaded (`.npz`).
5. **Outputs** — `metrics` turns trajectories into CIR / CMR; `TrajectoryStats` and
   `snapshot` implement the interval convention. `project()` wraps steps 1–5 and returns a
   tidy table.

## Design choices

* **Exact reproducibility.** Parameters are stored and read at full precision
  (`float_precision="round_trip"`), and every random draw comes from a generator seeded by
  (seed, coverage, simulation number), so a result never depends on what else is run.
  Regression tests (`pytest`) check regenerated ensembles and tables against
  `paper/reference`.
* **One design for every unit.** The national model and the 36 state/UT models share the
  equations, the sampling of `r_H` and `θ` (redrawn every year), the ensemble size (1000)
  and the seeds; they differ only in fitted parameters and demography.
* **Small, typed objects.** Parameters and inputs are dataclasses with JSON / CSV
  persistence, so a model is fully described by two files.
* **Paper code at arm's length.** `hpvcc.paper` only *uses* the core API. Deleting it would
  not affect the model.
* **One drawing code for paper and users.** Every figure type of the paper is a function
  of data in `hpvcc.viz.figures`; the paper subpackage passes it the paper's simulations,
  `Projection.plot_bars()` / `plot_fit()` pass it yours. Your figures therefore look exactly
  like the paper's, and with the paper's settings they *are* the paper's figures (checked
  pixel by pixel). Figures are saved inside their style context, because matplotlib reads
  the savefig settings (tight bounding box, padding) at save time.
* **Caching by content.** Paper experiments are cached in `results/paper/simulations/` and
  re-run automatically when the parameters or ensemble size change.

## Extending the model

**A new scenario type** (e.g. coverage that ramps up): add the time dependence to
`rhs_screening` (the coverage is read as `a = 0 if t < screening_start else coverage`),
or pass a callable coverage and evaluate it there; `Scenario` carries the settings.

**Other uncertainty assumptions for the screened pathway:** subclass `DrawScheme`:

```python
from hpvcc.simulation import DrawScheme
import numpy as np

class NarrowDraws(DrawScheme):
    name = "narrow"
    def draws(self, params, scenario_index, sim_index, n_steps):
        rng = np.random.default_rng(1000 * scenario_index + sim_index)
        return np.full(n_steps, rng.uniform(0.8, 1.0)), np.full(n_steps, rng.uniform(0, params.gamma / 20))
    def describe(self):
        return {"scheme": self.name}
```

**A new compartment or flow:** extend `COMPARTMENTS` in `constants.py`, both right-hand
sides in `model.py` (and `_FIELDS` if a parameter is added to `ModelParameters`), the
initial state built in `data/states.py` / `data/national.py`, then run
`hpvcc prepare-inputs`. The column indices of `C` (6) and `C_D` (9) are used by the
ensemble code and metrics; keep them or update `IDX` users. Results will no longer match
`paper/reference`, so the slow regression tests will (correctly) fail.

**A new unit:** add it to `REGIONS` in `units.py`, add its population rows to the
projection workbook and its prevalence to `data/raw/state_initial_hpv_prevalence.csv`,
add a parameter row to `params/states_final.csv`, and run `hpvcc prepare-inputs`.

## Tests

```bash
pytest -m "not slow"     # inputs, model, metrics, API, CLI, uncertainty band (~1 min)
pytest                   # + regenerated ensembles and LaTeX tables vs paper/reference
```
