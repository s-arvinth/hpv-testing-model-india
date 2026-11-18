# Using the model

Everything below runs from the fitted parameters and inputs shipped in the repository; no
other data are needed. Install first (`pip install -e .`, see the [README](../README.md#installation)).

- [Command line](#command-line)
- [Python: projections](#python-projections)
- [Python: the model directly](#python-the-model-directly)
- [Stochastic ensembles](#stochastic-ensembles)
- [Parameter uncertainty](#parameter-uncertainty)
- [Model fit](#model-fit)
- [Plots](#plots)
- [Performance](#performance)

## Command line

```text
hpvcc units                         list the national and 36 state/UT models
hpvcc params UNIT                   fitted parameters with symbols and meanings
hpvcc fit UNIT [--plot FILE]        fitted model vs observed rates, 2012–2019
hpvcc run UNIT [options]            project CIR and CMR, 2003–2070
hpvcc prepare-inputs                rebuild data/processed from data/raw
hpvcc paper ...                     reproduce the paper (see paper/README.md)
```

Unit names are case-insensitive, and `_` or `&` may be used (`tamil_nadu`,
`"Andaman & Nicobar"`).

### `hpvcc run`

| Option | Meaning | Default |
|---|---|---|
| `--coverage C [C ...]` | screening coverage(s) in [0, 1]; 0 = no screening | `0` |
| `--test T [T ...]` | `hpv-dna`, `lbc`, `via`, `pap` | `hpv-dna` |
| `--sensitivity H [H ...]` | test sensitivity H⁺ instead of `--test` | |
| `--start YEAR` | first year of screening | `2026` |
| `--sims N` | stochastic simulations per screened scenario | 1000 |
| `--deterministic` | single run at the midpoints of the r_H and θ ranges, no interval | off |
| `--seed S` | seed of the random draws | `0` |
| `--set NAME=VALUE ...` | change parameters | |
| `--params-file FILE` | parameters from your own CSV | |
| `--years Y [Y ...]` | years in the printed summary | 2030 … 2070 |
| `--out FILE.csv` | write the annual table | |
| `--plot FILE.png` | save CIR and CMR over time (`.png`, `.pdf`, `.svg`) | |
| `--bars FILE.png` | save the paper's coverage figure (one test) | |
| `--bar-years Y [Y ...]` | years of `--bars` | 2040 2050 2060 |
| `--workers N` (before `run`) | parallel processes | all cores − 2 |

Every coverage is combined with every test. Examples:

```bash
hpvcc run India --coverage 0 0.7                             # the paper's national numbers
hpvcc run India --coverage 0 0.25 0.5 0.7 0.9 1 --bars india.png   # = Fig 3 of the paper
hpvcc run Kerala --coverage 0 0.7                            # a state: the paper's numbers too
hpvcc run India --coverage 0.75 --test pap via lbc hpv-dna --sims 200
hpvcc run "Uttar Pradesh" --coverage 0 0.5 0.9 --start 2030 --plot up.png
hpvcc run Bihar --coverage 0 0.7 --set beta=0.5 r_I=0.6 --deterministic
hpvcc --workers 8 run Assam --coverage 0.7 --out assam.csv
```

With the defaults (1000 simulations, seed 0), runs reproduce the paper's numbers exactly,
for India and for every state/UT: the random draws of a scenario depend only on the seed,
its coverage (on the 5 % grid 0, 0.05, …, 1) and the simulation number, not on which other
scenarios are run with it.

## Python: projections

```python
from hpvcc import project

proj = project("India", coverage=[0, 0.5, 0.7], test="hpv-dna")
```

`project()` accepts the same choices as the command line: `coverage`, `test` or `H_pos`,
`n_sims`, `deterministic`, `params`, `params_file`, `changes`, `screening_start`, `seed`,
`workers`. It returns a [`Projection`](../src/hpvcc/projection.py):

| Attribute / method | Returns |
|---|---|
| `proj.table` | DataFrame, one row per scenario and year: `CIR`, `CIR_lo95`, `CIR_hi95`, `CMR`, …, `Cases`, `Deaths`, `Women` |
| `proj.summary(years)` | mean CIR / CMR at the given years, one row per scenario |
| `proj.at(2050)` | all scenarios in one year |
| `proj.elimination_year()` | first year from the start of screening with mean CIR < 4 |
| `proj.plot(path)` | CIR and CMR over time with 95 % bands |
| `proj.plot_bars(path, years)` | CIR and CMR by coverage at snapshot years (the paper's coverage figure) |
| `proj.snapshot(coverage, year)` | CIR, CMR and 95 % half-widths in one year, as in the figures |
| `proj.to_csv(path)` | the annual table |
| `proj.model`, `proj.ensemble` | the model used and the raw ensemble (if any) |

Rates are per 100,000 women; `Cases`, `Deaths` and `Women` are numbers of women.

## Python: the model directly

```python
from hpvcc import HPVCervicalCancerModel
from hpvcc.metrics import rates

model = HPVCervicalCancerModel.for_unit("India")
traj = model.simulate()                    # (69 years, 10 compartments), no screening
cir, cmr = rates(traj, model.params, model.inputs.denominator)

model.params                               # ModelParameters (frozen dataclass)
model.inputs.population["S30"]             # projected women aged 30–65
model.inputs.inflow(2030.5)                # recruitment Λ(t)
model.inputs.y0                            # initial state, 2002

what_if = model.with_params(gamma=2 * model.params.gamma).with_screening_start(2030)

# other initial conditions: build the model from parameters + modified inputs
more_infected = model.inputs.with_y0(I=1.2 * model.inputs.y0[5])
custom = HPVCervicalCancerModel(model.params, more_infected)
```

A screened scenario needs the two sampled rates of the screened pathway, which no
parameter set stores: `model.simulate(0.7, r_H, theta)`, e.g.
`model.simulate(0.7, *model.midpoint_rates())`; or run an ensemble (below).

## Stochastic ensembles

```python
from hpvcc.simulation import Scenario, StochasticEnsemble, YearlyDraws, coverage_sweep

scenarios = coverage_sweep([0.25, 0.5, 0.75], H_pos=0.95) + [Scenario.for_test(0.75, "via")]
ens = StochasticEnsemble(model, scenarios, n_sims=300, draws=YearlyDraws(seed=7),
                         keep_full=True).run()
ens.snapshot(ens.index(coverage=0.5), 2050)       # CIR, CMR and half-widths in 2050
ens.stats(0, "H_T")                               # mean / 2.5 / 97.5 % trajectories of any compartment
ens.save("my_ensemble.npz")
```

`YearlyDraws` redraws `r_H ~ U(r_I, 1)` and `θ ~ U(0, γ)` every year; to try another assumption, subclass
[`DrawScheme`](../src/hpvcc/simulation.py) and implement `draws()` and `describe()`.

## Parameter uncertainty

```python
from hpvcc.uncertainty import ParameterSetEnsemble

res = ParameterSetEnsemble.from_csv(HPVCervicalCancerModel.for_unit("India")).run()
res.to_frame()          # year, mean, 2.5 / 97.5 percentiles and best-set CIR and CMR
```

Any table of parameter sets works (`ParameterSetEnsemble(model, my_dataframe)`), e.g. to
propagate your own priors through the no-screening model.

## Model fit

```python
from hpvcc.observed import ObservedRates, model_fit
from hpvcc.viz.plots import plot_fit

model_fit("India")                           # Year, Model_CIR, Observed_CIR, Model_CMR, Observed_CMR
ObservedRates.load("Kerala").cir             # observed CIR 2012–2019 of a state
plot_fit("India", "fit.png")                 # Fig 1 (a state: its panel of Fig 2)
plot_fit("India", "fit2.png", model=model.with_params(gamma=2 * model.params.gamma))
```

The observed rates are the ICMR-NCDIR / NINE registry estimates in
[`data/registry/`](../data/registry/); they are shown, never fitted, by this package.

## Plots

Every figure the package draws uses the paper's figure code, so a figure of your own run
looks exactly like the corresponding figure of the paper — same layout, fonts (Fira Sans),
sizes, colours (each state's region hue, viridis for India, the cool palette for the four
tests), legends and ticks:

| Function | Figure type | As in the paper |
|---|---|---|
| `proj.plot()` / `plot_projection` | CIR and CMR over time, 95 % bands | time-series panel of Fig S3 |
| `proj.plot_bars()` / `plot_coverage_bars` | CIR and CMR by coverage at three years | Fig 3, Figs S11–S16 |
| `plot_fit` | fitted model vs observed, 2012–2019 | Fig 1, Fig 2, Figs S5–S10 |
| `plot_compartments` | every compartment over time | time-series panel of Fig S3 |
| `hpvcc.viz.figures.test_bars` | several tests: detected infections, CIR/CMR bars | Fig S3 |
| `hpvcc.viz.figures.uncertainty_bars` | parameter-uncertainty band | Fig S4 |

With the paper's settings the output is the paper's figure: e.g.
`project("India", coverage=[0, .25, .5, .7, .9, 1]).plot_bars("india.pdf")` is
pixel-identical to Fig 3, `plot_fit("India", ...)` to Fig 1, and the same calls for a
state to its panels of Figs S11–S16 and Fig 2. Functions save to the path given (format from the extension) or, with
no path, return the open figure.

**Background.** Figures are saved on a white background. To save them without one
(transparent figure and axes; legend boxes kept), as for slides or dark pages:

```python
proj.plot("india.png", transparent=True)                    # one figure

from hpvcc.viz.style import set_transparent, transparent_background
set_transparent(True)                                    # every following figure
with transparent_background():                           # figures inside the block
    proj.plot_bars("india_bars.png")
```

On the command line, `hpvcc --no-background run ... --plot india.png`; for the paper,
`hpvcc --no-background paper figure fig3` writes to `results/paper/figures_transparent/`
(and the model diagrams in their white-line version). `HPVCC_TRANSPARENT=1` sets the
default for a whole session.

`hpvcc.viz.style.figure_style("paper")` applies the paper's style to figures of your own:

```python
import matplotlib.pyplot as plt
from hpvcc.viz.style import figure_style, save

with figure_style("paper"):
    fig, ax = plt.subplots(figsize=(14, 5.2))
    ...
    save(fig, "my_figure.pdf")        # save inside the block: savefig settings are part of the style
```

## Performance

A deterministic run takes about 0.05 s. A stochastic simulation costs 68 one-year ODE
integrations (about 25 ms on one core): 1000 simulations of one scenario
take about 30 s on a single core and about a second on 40 cores. Use `--sims` / `n_sims` for quick looks
and `--workers` to set the number of processes.
