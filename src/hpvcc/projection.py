"""
High-level projections: one call from a unit and a screening scenario to a tidy table
of yearly crude incidence and mortality rates.

>>> from hpvcc import project
>>> proj = project("India", coverage=[0, 0.7])           # no screening vs 70 % HPV-DNA testing
>>> proj.summary()                                       # CIR / CMR at 2030, 2040, ..., 2070
>>> proj.table                                           # every year, with 95 % intervals
>>> proj.plot("india.png")                               # rates over time
>>> proj.plot_bars("india_bars.png")                     # paper-style bars at 2040/2050/2060

How screened scenarios are simulated
------------------------------------
By default each screened scenario is a stochastic ensemble with the draw scheme and size
used in the paper (:func:`hpvcc.simulation.default_draws`, 1000 simulations): the clearance
and progression rates of screen-detected infection, r_H and θ, are redrawn every year
(r_H ~ U(r_I, 1), θ ~ U(0, γ)) and the 95 % interval
is mean ± half the spread of the 2.5/97.5-percentile trajectories. With
``deterministic=True`` a single run uses the midpoints of those ranges instead (fast, no
interval). The no-screening scenario is deterministic either way.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd

from .constants import PER_100K, SCREENING_START, SCREENING_TESTS, WHO_THRESHOLD
from .metrics import RateSnapshot, TrajectoryStats, death_counts, incidence_counts, snapshot
from .model import HPVCervicalCancerModel
from .parameters import ModelParameters, load_parameters
from .simulation import (EnsembleResult, Scenario, StochasticEnsemble, default_draws, default_n_sims,
                         sensitivity_of)
from .units import INDIA, resolve_unit

SUMMARY_YEARS = (2030, 2040, 2050, 2060, 2070)
_TEST_NAME = {round(h, 2): name for name, h in SCREENING_TESTS.items()}


def _as_list(x) -> list:
    if x is None:
        return []
    if isinstance(x, (str, bytes)) or not isinstance(x, Iterable):
        return [x]
    return list(x)


def _rate_rows(unit: str, scenario: Scenario, years: np.ndarray, denominator: np.ndarray,
               params: ModelParameters, C: TrajectoryStats, D: TrajectoryStats, method: str) -> pd.DataFrame:
    """Yearly rates of one scenario from trajectory statistics (the paper's interval convention)."""
    k = PER_100K / np.maximum(denominator[1:], 1e-12)
    im, il, ih = (incidence_counts(x, params) for x in (C.mean, C.low, C.high))
    dm, dl, dh = (death_counts(x) for x in (D.mean, D.low, D.high))
    hwi, hwd = np.abs(ih - il) / 2, np.abs(dh - dl) / 2
    return pd.DataFrame({
        "Unit": unit, "Coverage": scenario.coverage, "H_pos": scenario.H_pos,
        "Test": "No screening" if scenario.coverage == 0 else
                _TEST_NAME.get(round(scenario.H_pos, 2), f"H+={scenario.H_pos:g}"),
        "Year": years[1:].astype(int),
        "CIR": im * k, "CIR_lo95": np.maximum(im - hwi, 0) * k, "CIR_hi95": (im + hwi) * k,
        "CMR": dm * k, "CMR_lo95": np.maximum(dm - hwd, 0) * k, "CMR_hi95": (dm + hwd) * k,
        "Cases": im * 1000, "Deaths": dm * 1000,
        "Women": denominator[1:] * 1000, "Method": method,
    })


@dataclass
class Projection:
    """Result of :func:`project`.

    Attributes
    ----------
    unit : str
    scenarios : list of Scenario
    table : DataFrame
        One row per scenario and year (2003–2070): ``CIR``, ``CMR`` (per 100,000 women) with
        95 % intervals (``_lo95`` / ``_hi95``), ``Cases``, ``Deaths`` and ``Women`` (numbers).
    model : HPVCervicalCancerModel
        The model that produced it (with any parameter changes applied).
    ensemble : EnsembleResult or None
        The raw stochastic ensemble, if one was run.
    """

    unit: str
    scenarios: list[Scenario]
    table: pd.DataFrame
    model: HPVCervicalCancerModel
    ensemble: EnsembleResult | None = None
    settings: dict = field(default_factory=dict)
    stats: dict = field(default_factory=dict, repr=False)   #: (coverage, H_pos) -> (C stats, C_D stats)

    def at(self, year: int) -> pd.DataFrame:
        """All scenarios in one year."""
        return self.table[self.table.Year == int(year)].reset_index(drop=True)

    def summary(self, years: Sequence[int] = SUMMARY_YEARS, metrics=("CIR", "CMR")) -> pd.DataFrame:
        """Mean rates at a few years, one row per scenario (``CIR_2050``, ``CMR_2050``, ...)."""
        t = self.table[self.table.Year.isin(years)]
        wide = t.pivot_table(index=["Coverage", "Test"], columns="Year", values=list(metrics), sort=False)
        wide.columns = [f"{m}_{y}" for m, y in wide.columns]
        return wide.reset_index()

    def elimination_year(self, threshold: float = WHO_THRESHOLD) -> pd.DataFrame:
        """First year (from the start of screening) in which the mean CIR is below ``threshold``."""
        start = self.model.screening_start
        rows = []
        for (cov, test), g in self.table.groupby(["Coverage", "Test"], sort=False):
            below = g[(g.Year >= start) & (g.CIR < threshold)]
            rows.append({"Coverage": cov, "Test": test,
                         "First_year_CIR_below": int(below.Year.iloc[0]) if len(below) else None})
        out = pd.DataFrame(rows)
        out["First_year_CIR_below"] = out["First_year_CIR_below"].astype("Int64")
        return out

    def to_csv(self, path: str | Path) -> Path:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        self.table.to_csv(path, index=False)
        return path

    def snapshot(self, coverage: float, year: int, H_pos: float | None = None) -> RateSnapshot:
        """CIR / CMR and their 95 % half-widths in one year, computed exactly as in the
        paper's figures (no screening: computed if it was not projected)."""
        key = self._key(coverage, H_pos)
        if key is None and coverage == 0:
            traj = self.model.simulate()
            C = TrajectoryStats(traj[:, 6], traj[:, 6], traj[:, 6])
            D = TrajectoryStats(traj[:, 9], traj[:, 9], traj[:, 9])
        elif key is None:
            raise KeyError(f"no scenario with coverage={coverage}, H_pos={H_pos}")
        else:
            C, D = self.stats[key]
        return snapshot(C, D, self.model.years, int(year), self.model.params, self.model.inputs.denominator)

    def _key(self, coverage, H_pos):
        for c, h in self.stats:
            if round(c, 2) == round(coverage, 2) and (coverage == 0 or H_pos is None or round(h, 2) == round(H_pos, 2)):
                return (c, h)
        return None

    def coverage_snapshots(self, years):
        """``({coverage: {year: (CIR, hw, CMR, hw)}}, coverages)`` for the coverage figure."""
        tests = {round(s.H_pos, 2) for s in self.scenarios if s.coverage > 0}
        if len(tests) > 1:
            raise ValueError("the coverage figure needs a single screening test; "
                             "project one test or use plot() to compare tests")
        covs = sorted({0.0, *(s.coverage for s in self.scenarios)})
        return {a: {Y: self.snapshot(a, Y).as_tuple() for Y in years} for a in covs}, covs

    def plot(self, path: str | Path | None = None, **kwargs):
        """CIR and CMR over time, one line (and 95 % band) per scenario. See
        :func:`hpvcc.viz.plots.plot_projection`."""
        from .viz.plots import plot_projection
        return plot_projection(self, path, **kwargs)

    def plot_bars(self, path: str | Path | None = None, years=(2040, 2050, 2060), **kwargs):
        """CIR and CMR by coverage at snapshot years, as in the paper's coverage figures
        (Fig 3, Figs S11–S16). See :func:`hpvcc.viz.plots.plot_coverage_bars`."""
        from .viz.plots import plot_coverage_bars
        return plot_coverage_bars(self, path, years, **kwargs)


def project(unit: str, coverage: float | Sequence[float] = 0.0, test: str | Sequence[str] | None = None,
            H_pos: float | Sequence[float] | None = None, *, n_sims: int | None = None,
            deterministic: bool = False, params: ModelParameters | None = None,
            params_file: str | Path | None = None, changes: dict | None = None,
            screening_start: float = SCREENING_START, seed: int = 0,
            workers: int | None = None) -> Projection:
    """Project CIR and CMR for one unit under one or more screening scenarios.

    Parameters
    ----------
    unit : str
        ``"India"`` or a state/UT (``hpvcc units`` lists them; case and ``_`` are ignored).
    coverage : float or list of float
        Screening coverage α in [0, 1] (0 = no screening). Every coverage is combined with
        every test.
    test : str or list of str, optional
        Screening test(s) by name: ``"hpv-dna"`` (default), ``"lbc"``, ``"via"``, ``"pap"``.
    H_pos : float or list of float, optional
        Test sensitivity given directly (instead of ``test``).
    n_sims : int, optional
        Stochastic simulations per screened scenario (default 1000, as in the paper).
    deterministic : bool
        Use a single run at the midpoints of the r_H and θ ranges instead of an ensemble.
    params : ModelParameters, optional
        A complete parameter set to use instead of the fitted one.
    params_file : path, optional
        A CSV of parameters (see :func:`hpvcc.parameters.load_parameters`).
    changes : dict, optional
        Individual parameter changes applied last, e.g. ``{"beta": 0.4}``.
    screening_start : float
        First year of screening (default 2026).
    seed : int
        Seed of the random draws (default 0, as in the paper).
    workers : int, optional
        Parallel processes for ensembles (default: all cores but two, at most 48).
    """
    unit = resolve_unit(unit)
    if params is None:
        params = load_parameters(unit, file=params_file) if params_file else load_parameters(unit)
    model = HPVCervicalCancerModel.for_unit(unit, params=params, screening_start=screening_start)
    if changes:
        model = model.with_params(**changes)

    sens = [float(h) for h in _as_list(H_pos)] or [sensitivity_of(t) for t in _as_list(test)] \
        or [model.params.H_pos]
    covs = [float(c) for c in _as_list(coverage)] or [0.0]
    for c in covs:
        if not 0.0 <= c <= 1.0:
            raise ValueError(f"coverage must be within [0, 1], got {c}")
    scenarios: list[Scenario] = []
    for c in covs:
        for h in ([sens[0]] if c == 0.0 else sens):    # one no-screening arm is enough
            scenarios.append(Scenario(c, h))

    p = model.params
    years, N = model.years, model.inputs.denominator
    frames, ensemble, stats = [], None, {}
    screened = [s for s in scenarios if s.coverage > 0]

    if screened and not deterministic:
        n = int(n_sims or default_n_sims(unit))
        ensemble = StochasticEnsemble(model, screened, n, default_draws(unit, seed, screened), workers=workers).run()

    for s in scenarios:
        if s.coverage == 0.0:
            traj = model.simulate()
            C = TrajectoryStats(traj[:, 6], traj[:, 6], traj[:, 6])
            D = TrajectoryStats(traj[:, 9], traj[:, 9], traj[:, 9])
            frames.append(_rate_rows(unit, s, years, N, p, C, D, "deterministic"))
            stats[(s.coverage, s.H_pos)] = (C, D)
        elif deterministic:
            m = model.with_params(H_pos=s.H_pos)
            traj = m.simulate(s.coverage, *m.midpoint_rates())
            C = TrajectoryStats(traj[:, 6], traj[:, 6], traj[:, 6])
            D = TrajectoryStats(traj[:, 9], traj[:, 9], traj[:, 9])
            frames.append(_rate_rows(unit, s, years, N, p, C, D, "deterministic (mid-range r_H, θ)"))
            stats[(s.coverage, s.H_pos)] = (C, D)
        else:
            i = ensemble.index(coverage=s.coverage, H_pos=s.H_pos)
            C, D = ensemble.stats(i, "C"), ensemble.stats(i, "C_D")
            frames.append(_rate_rows(unit, s, years, N, p, C, D, f"stochastic ({ensemble.n_sims} simulations)"))
            stats[(s.coverage, s.H_pos)] = (C, D)

    settings = {"n_sims": None if ensemble is None else ensemble.n_sims, "deterministic": deterministic,
                "screening_start": screening_start, "changes": dict(changes or {}),
                "params_file": None if params_file is None else str(params_file)}
    return Projection(unit, scenarios, pd.concat(frames, ignore_index=True), model, ensemble, settings, stats)
