"""Refresh paper/reference from results/paper after an INTENDED change of results.

    hpvcc paper all && python paper/update_reference.py

Copies the LaTeX tables and figures, and rewrites reference_numbers.csv (CIR, CMR and
95 % half-widths of every unit at coverage 0–100 % in 5 % steps, every fifth year
2030–2070) from the cached coverage sweeps. Review the diff before committing.
"""
import shutil

import pandas as pd

from hpvcc.paper.experiments import Experiment
from hpvcc.paths import FIGURES_DIR, REFERENCE_DIR, TABLES_DIR
from hpvcc.units import INDIA, STATES

for name in ("paper_tables", "paper_figures"):
    shutil.rmtree(REFERENCE_DIR / name, ignore_errors=True)
shutil.copytree(TABLES_DIR / "latex", REFERENCE_DIR / "paper_tables")
shutil.copytree(FIGURES_DIR, REFERENCE_DIR / "paper_figures")

rows = []
for unit in [*STATES, INDIA]:
    res = Experiment.get("india_coverage_sweep" if unit == INDIA else "states_coverage_sweep").result(unit)
    for i, sc in enumerate(res.scenarios):
        for year in range(2030, 2071, 5):
            s = res.snapshot(i, year)
            rows.append((unit, year, int(round(sc.coverage * 100)), round(s.cir, 3), round(s.cmr, 3),
                         round(s.cir_halfwidth, 3), round(s.cmr_halfwidth, 3)))
pd.DataFrame(rows, columns=["State", "Year", "Coverage_pct", "CIR", "CMR", "CIR_95CI_hw", "CMR_95CI_hw"]).to_csv(
    REFERENCE_DIR / "reference_numbers.csv", index=False)
print("reference updated:", REFERENCE_DIR)
