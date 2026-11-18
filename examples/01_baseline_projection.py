"""
01 — No-screening projection for India.

Loads the fitted parameters and inputs of the national model, integrates it from 2002 to
2070 and prints the crude incidence (CIR) and mortality (CMR) rates per 100,000 women;
then plots every compartment of the trajectory. Any state/UT works the same way
(`project("Kerala")`, `hpvcc units` lists them).

    python examples/01_baseline_projection.py
"""
from pathlib import Path

from hpvcc import HPVCervicalCancerModel, project
from hpvcc.viz.plots import plot_compartments

OUT = Path("results/examples")

proj = project("India")                                   # coverage 0 = no screening
print("India: no-screening projection")
print(proj.summary().round(2).to_string(index=False))

# Every compartment of the national trajectory (thousands of women)
model = HPVCervicalCancerModel.for_unit("India")
print(f"\nfigure -> {plot_compartments(model, path=OUT / '01_india_compartments.png')}")
