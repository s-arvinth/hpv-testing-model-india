"""
Coverage-sweep snapshot bars (Fig 3 for India, Figs S11–S16 for the states/UTs), drawn by
:func:`hpvcc.viz.figures.coverage_bars` from the paper's coverage sweeps.
"""
from __future__ import annotations

from pathlib import Path

from ...units import INDIA
from ...viz.figures import SWEEP_COVERAGES, SWEEP_YEARS, coverage_bars, viridis_sweep_cmap  # noqa: F401
from ..experiments import Experiment


def coverage_sweep_bars(unit: str, out: Path, years=SWEEP_YEARS, coverages=SWEEP_COVERAGES,
                        experiment: str | None = None) -> Path:
    exp = Experiment.get(experiment or ("india_coverage_sweep" if unit == INDIA else "states_coverage_sweep"))
    res = exp.result(unit)
    snap = {a: {Y: res.snapshot(res.index(coverage=a), Y).as_tuple() for Y in years} for a in coverages}
    return coverage_bars(unit, snap, out, years, coverages)
