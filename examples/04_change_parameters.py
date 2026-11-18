"""
04 — What-if: run the model with changed parameters.

Three ways to use other parameter values:
  * individual changes:       project(..., changes={"beta": 0.4})
  * a full parameter object:  project(..., params=load_parameters("India").replace(...))
  * a CSV of your own:        project(..., params_file="my_params.csv")
The same options exist on the command line: --set beta=0.4, --params-file my_params.csv.

    python examples/04_change_parameters.py
"""
from pathlib import Path

import pandas as pd

from hpvcc import HPVCervicalCancerModel, load_parameters, project
from hpvcc.viz.plots import plot_fit

OUT = Path("results/examples")
OUT.mkdir(parents=True, exist_ok=True)
unit = "India"
p = load_parameters(unit)
print(f"{unit}: fitted beta = {p.beta:.3f}, gamma = {p.gamma:.5f}, r_I = {p.r_I:.3f}")

runs = {
    "fitted": project(unit, coverage=[0, 0.7], deterministic=True),
    "transmission +20 %": project(unit, coverage=[0, 0.7], deterministic=True, changes={"beta": 1.2 * p.beta}),
    "progression x2": project(unit, coverage=[0, 0.7], deterministic=True, changes={"gamma": 2 * p.gamma}),
    "screening from 2030": project(unit, coverage=[0, 0.7], deterministic=True, screening_start=2030),
}
rows = []
for name, proj in runs.items():
    s = proj.summary(years=(2050,))
    rows.append({"run": name, "CIR 2050, no screening": s.CIR_2050.iloc[0], "CIR 2050, 70 %": s.CIR_2050.iloc[1]})
print(pd.DataFrame(rows).round(2).to_string(index=False))

# How a change moves the fit to the observed 2012–2019 rates (layout of the paper's Fig 1)
changed = HPVCervicalCancerModel.for_unit(unit).with_params(gamma=2 * p.gamma)
print(f"\nfigure -> {plot_fit(unit, OUT / '04_india_fit_progression_x2.png', model=changed)}")

# Your own parameter file: same columns as params/india_final.csv (one row, or one per unit)
my_file = OUT / "04_my_parameters.csv"
pd.DataFrame([p.replace(gamma=1.2 * p.gamma).as_dict()]).to_csv(my_file, index=False)   # progression +20 %
mine = project(unit, coverage=[0, 0.7], deterministic=True, params_file=my_file)
print(f"\nwith {my_file}:")
print(mine.summary(years=(2030, 2050, 2070)).round(2).to_string(index=False))
