"""
05 — The building blocks underneath `project`.

Model = parameters + inputs; scenarios; a draw scheme for the uncertain rates of the
screened pathway; an ensemble runner; summaries of compartments and rates.

    python examples/05_low_level_api.py
"""
import numpy as np

from hpvcc.constants import COMPARTMENTS
from hpvcc.metrics import rates
from hpvcc.model import HPVCervicalCancerModel
from hpvcc.simulation import Scenario, StochasticEnsemble, YearlyDraws

model = HPVCervicalCancerModel.for_unit("India")
print("initial state (2002, thousands):", {c: round(float(v), 1) for c, v in zip(COMPARTMENTS, model.inputs.y0)})

# Deterministic no-screening trajectory: array (years, 10 compartments)
traj = model.simulate()
cir, cmr = rates(traj, model.params, model.inputs.denominator)
print(f"no screening: CIR 2050 = {cir[list(model.years[1:]).index(2050)]:.2f}")

# Stochastic ensemble: 70 % coverage with VIA and with HPV-DNA testing, 100 simulations each,
# r_H ~ U(r_I, 1) and theta ~ U(0, gamma) redrawn every year
scenarios = [Scenario.for_test(0.7, "via"), Scenario.for_test(0.7, "hpv-dna")]
ens = StochasticEnsemble(model, scenarios, n_sims=100, draws=YearlyDraws(seed=1),
                         keep_full=True).run()
for i, sc in enumerate(ens.scenarios):
    snap = ens.snapshot(i, 2050)
    H_T = ens.stats(i, "H_T")                     # any compartment, mean and 2.5/97.5 percentiles
    print(f"{sc.label}: CIR 2050 = {snap.cir:.2f} ± {snap.cir_halfwidth:.2f}; "
          f"screen-detected infections 2030 = {H_T.mean[list(model.years).index(2030)]:.1f} thousand")
