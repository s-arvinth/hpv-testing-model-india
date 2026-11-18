"""
Stochastic screening ensembles.

Two rates of the screened pathway are uncertain and are sampled in every simulation, anew
for every simulated year: the clearance rate of screen-detected infection,
``r_H ~ U(r_I, 1)``, and its progression rate to cancer, ``θ ~ U(0, γ)``. Everything else
is fixed at the fitted values. Because H_T is empty without screening, the no-screening
arm is deterministic and the ensemble spread only reflects the screened pathway.

Every unit — the national model and the 36 state/UT models — uses the same scheme
(:class:`YearlyDraws`) and the same ensemble size (:data:`DEFAULT_N_SIMS` = 1000).
Simulation ``j`` of a scenario draws from ``default_rng([seed, key, j])``, where ``key``
is the scenario's coverage on a 5 % grid: a result at a given coverage does not depend on
which other scenarios are run with it, and scenarios at the same coverage (e.g. different
tests) share their draws (common random numbers).
"""
from __future__ import annotations

import json
import os
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from multiprocessing import Pool
from pathlib import Path
from typing import Sequence

import numpy as np

from .constants import IDX
from .metrics import RateSnapshot, TrajectoryStats, snapshot
from .model import HPVCervicalCancerModel, pack
from .parameters import ModelParameters

# --------------------------------------------------------------------------------------
# Scenarios
# --------------------------------------------------------------------------------------


@dataclass(frozen=True)
class Scenario:
    """One screening scenario: programme coverage and test sensitivity."""

    coverage: float
    H_pos: float

    @property
    def label(self) -> str:
        return f"coverage={self.coverage:.2f}, H+={self.H_pos:.2f}"


    @classmethod
    def for_test(cls, coverage: float, test: str) -> "Scenario":
        """Scenario for a named screening test (``"VIA"``, ``"pap"``, ``"hpv-dna"``, ...)."""
        return cls(float(coverage), sensitivity_of(test))


def sensitivity_of(test: str) -> float:
    """Sensitivity H+ of a screening test given by name or alias."""
    from .constants import SCREENING_TESTS, TEST_ALIASES
    key = TEST_ALIASES.get(test.strip().lower(), test.strip())
    for name, h in SCREENING_TESTS.items():
        if name.lower() == key.lower():
            return h
    raise KeyError(f"unknown screening test {test!r}; choose from {', '.join(TEST_ALIASES)}")


def coverage_sweep(coverages, H_pos: float) -> list[Scenario]:
    return [Scenario(float(c), float(H_pos)) for c in coverages]


def sensitivity_sweep(sensitivities, coverage: float) -> list[Scenario]:
    return [Scenario(float(coverage), float(h)) for h in sensitivities]


# --------------------------------------------------------------------------------------
# Random-draw schemes
# --------------------------------------------------------------------------------------


class DrawScheme(ABC):
    """How (r_H, θ) are drawn for each simulation of each scenario."""

    name: str = ""

    def prepare(self, params: ModelParameters, n_scenarios: int, n_sims: int, n_steps: int) -> None:
        """Hook for schemes that must generate all draws up front (in loop order)."""

    @abstractmethod
    def draws(self, params: ModelParameters, scenario_index: int, sim_index: int,
              n_steps: int) -> tuple[np.ndarray, np.ndarray]:
        """``(r_H, θ)`` arrays of length ``n_steps`` for one simulation."""

    @abstractmethod
    def describe(self) -> dict:
        """JSON-serialisable description (stored with the results)."""


class YearlyDraws(DrawScheme):
    """New ``r_H ~ U(r_I, 1)`` and ``θ ~ U(0, γ)`` every year.

    Simulation ``j`` of scenario ``i`` uses ``numpy.random.default_rng([seed, key_i, j])``;
    ``key_i`` is ``scenario_keys[i]`` (see :func:`coverage_keys`) or, without keys, ``i``.
    """

    name = "yearly"

    def __init__(self, seed: int = 0, scenario_keys: Sequence[int] | None = None):
        self.seed = int(seed)
        self.keys = None if scenario_keys is None else [int(k) for k in scenario_keys]

    def seed_of(self, scenario_index: int, sim_index: int) -> list[int]:
        key = scenario_index if self.keys is None else self.keys[scenario_index]
        return [self.seed, key, sim_index]

    def draws(self, params, scenario_index, sim_index, n_steps):
        rng = np.random.default_rng(self.seed_of(scenario_index, sim_index))
        r_H = np.empty(n_steps)
        theta = np.empty(n_steps)
        for k in range(n_steps):
            r_H[k] = rng.uniform(params.r_I, 1.0)
            theta[k] = rng.uniform(0.0, params.gamma)
        return r_H, theta

    def describe(self):
        d = {"scheme": self.name, "seed": self.seed, "rng": "default_rng([seed, scenario_key, sim])",
             "r_H": "U(r_I, 1) every year", "theta": "U(0, gamma) every year"}
        if self.keys is not None:
            d["scenario_keys"] = self.keys
        return d


def coverage_keys(scenarios: Sequence["Scenario"], grid: float = 0.05) -> list[int]:
    """Seed key of each scenario: its coverage on a 5 % grid (0 → 0, 0.05 → 1, …, 1 → 20);
    coverages off the grid get keys from 1000 upwards."""
    keys = []
    for i, s in enumerate(scenarios):
        k = s.coverage / grid
        keys.append(int(round(k)) if abs(k - round(k)) < 1e-9 else 1000 + i)
    return keys


DRAW_SCHEMES = {cls.name: cls for cls in (YearlyDraws,)}

#: Simulations per screened scenario, for every unit.
DEFAULT_N_SIMS = 1000
DEFAULT_SEED = 0


def default_draws(unit: str | None = None, seed: int = DEFAULT_SEED,
                  scenarios: Sequence["Scenario"] | None = None) -> DrawScheme:
    """The draw scheme of the paper (the same for every unit): :class:`YearlyDraws`, keyed
    by coverage when the ``scenarios`` are given."""
    return YearlyDraws(seed, None if scenarios is None else coverage_keys(scenarios))


def default_n_sims(unit: str | None = None) -> int:
    return DEFAULT_N_SIMS


def make_draw_scheme(spec: dict, scenarios: Sequence["Scenario"] | None = None) -> DrawScheme:
    """Build a scheme from a config entry such as ``{"scheme": "yearly", "seed": 0}``."""
    spec = dict(spec)
    cls = DRAW_SCHEMES[spec.pop("scheme")]
    if cls is YearlyDraws and scenarios is not None:
        spec.setdefault("scenario_keys", coverage_keys(scenarios))
    return cls(**spec)


# --------------------------------------------------------------------------------------
# Results container
# --------------------------------------------------------------------------------------


@dataclass
class EnsembleResult:
    """Output of a stochastic ensemble for one unit.

    ``C`` and ``C_D`` have shape ``(n_scenarios, n_sims, n_years)``; ``full`` (optional)
    holds every compartment, shape ``(n_scenarios, n_sims, n_years, 10)``.
    """

    unit: str
    scenarios: list[Scenario]
    years: np.ndarray
    denominator: np.ndarray
    params: ModelParameters
    C: np.ndarray
    C_D: np.ndarray
    full: np.ndarray | None = None
    meta: dict = field(default_factory=dict)

    # ------------------------------------------------------------------ access
    @property
    def n_sims(self) -> int:
        return self.C.shape[1]

    def index(self, *, coverage: float | None = None, H_pos: float | None = None) -> int:
        """Scenario index by coverage and/or sensitivity (rounded to 2 dp)."""
        for i, s in enumerate(self.scenarios):
            if ((coverage is None or round(s.coverage, 2) == round(coverage, 2)) and
                    (H_pos is None or round(s.H_pos, 2) == round(H_pos, 2))):
                return i
        raise KeyError(f"no scenario with coverage={coverage}, H_pos={H_pos}")

    def stats(self, scenario_index: int, compartment: str = "C") -> TrajectoryStats:
        if compartment == "C":
            sims = self.C[scenario_index]
        elif compartment == "C_D":
            sims = self.C_D[scenario_index]
        else:
            if self.full is None:
                raise ValueError("this ensemble stored only C and C_D")
            sims = self.full[scenario_index, :, :, IDX[compartment]]
        return TrajectoryStats.of(sims, axis=0)

    def snapshot(self, scenario_index: int, year: int) -> RateSnapshot:
        """CIR/CMR (mean ± figure half-width) in ``year`` for one scenario."""
        return snapshot(self.stats(scenario_index, "C"), self.stats(scenario_index, "C_D"),
                        self.years, year, self.params, self.denominator)

    def subset(self, indices) -> "EnsembleResult":
        indices = list(indices)
        return EnsembleResult(self.unit, [self.scenarios[i] for i in indices], self.years,
                              self.denominator, self.params, self.C[indices], self.C_D[indices],
                              None if self.full is None else self.full[indices], dict(self.meta))

    # ------------------------------------------------------------------ persistence
    def save(self, path: Path) -> Path:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        arrays = dict(C=self.C, C_D=self.C_D, years=self.years, denominator=self.denominator,
                      coverage=np.array([s.coverage for s in self.scenarios]),
                      H_pos=np.array([s.H_pos for s in self.scenarios]))
        if self.full is not None:
            arrays["full"] = self.full
        meta = dict(self.meta, unit=self.unit, params=self.params.as_dict())
        np.savez_compressed(path, meta=np.array(json.dumps(meta)), **arrays)
        return path

    @classmethod
    def load(cls, path: Path) -> "EnsembleResult":
        z = np.load(path, allow_pickle=False)
        meta = json.loads(str(z["meta"]))
        scen = [Scenario(float(c), float(h)) for c, h in zip(z["coverage"], z["H_pos"])]
        return cls(meta["unit"], scen, z["years"], z["denominator"],
                   ModelParameters(**meta["params"]), z["C"], z["C_D"],
                   z["full"] if "full" in z.files else None, meta)


# --------------------------------------------------------------------------------------
# Ensemble runner
# --------------------------------------------------------------------------------------

_WORKER: dict = {}


def _init_worker(model: HPVCervicalCancerModel, scheme: DrawScheme, keep_full: bool):
    _WORKER.update(model=model, scheme=scheme, keep_full=keep_full)


def _run_one(task):
    scen_index, sim_index, scenario = task
    model: HPVCervicalCancerModel = _WORKER["model"]
    scheme: DrawScheme = _WORKER["scheme"]
    params = model.params
    ys = model.years
    nT = len(ys)
    r_H, theta = scheme.draws(params, scen_index, sim_index, nT - 1)
    P = pack(params.replace(H_pos=scenario.H_pos))
    y = np.array(model.inputs.y0, dtype=float)
    out = np.empty((nT, 10))
    out[0] = y
    for i in range(1, nT):
        y = model.step_screened(y, ys[i - 1], ys[i], scenario.coverage, r_H[i - 1], theta[i - 1], P)
        out[i] = y
    return out if _WORKER["keep_full"] else out[:, [6, 9]]


def default_workers() -> int:
    env = os.environ.get("HPVCC_WORKERS")
    if env:
        return max(1, int(env))
    return max(1, min(48, (os.cpu_count() or 2) - 2))


class StochasticEnsemble:
    """Run ``n_sims`` stochastic screening simulations for each scenario.

    Parameters
    ----------
    model : HPVCervicalCancerModel
    scenarios : list of Scenario
    n_sims : int
    draws : DrawScheme
    keep_full : bool
        Keep all 10 compartments (needed for compartment figures); otherwise only C and C_D.
    workers : int, optional
        Parallel processes (default: ``HPVCC_WORKERS`` or ``cpu_count - 2``, max 48).
    """

    def __init__(self, model: HPVCervicalCancerModel, scenarios: list[Scenario], n_sims: int,
                 draws: DrawScheme, keep_full: bool = False, workers: int | None = None):
        self.model = model
        self.scenarios = list(scenarios)
        self.n_sims = int(n_sims)
        self.draws = draws
        self.keep_full = keep_full
        self.workers = workers or default_workers()

    def run(self) -> EnsembleResult:
        nT = len(self.model.years)
        self.draws.prepare(self.model.params, len(self.scenarios), self.n_sims, nT - 1)
        tasks = [(si, j, sc) for si, sc in enumerate(self.scenarios) for j in range(self.n_sims)]
        args = (self.model, self.draws, self.keep_full)
        if self.workers > 1:
            with Pool(self.workers, initializer=_init_worker, initargs=args) as pool:
                outs = pool.map(_run_one, tasks, chunksize=max(1, len(tasks) // (self.workers * 8)))
        else:
            _init_worker(*args)
            outs = [_run_one(t) for t in tasks]
        arr = np.asarray(outs).reshape(len(self.scenarios), self.n_sims, nT, -1)
        if self.keep_full:
            full, C, CD = arr, arr[..., 6], arr[..., 9]
        else:
            full, C, CD = None, arr[..., 0], arr[..., 1]
        meta = {"n_sims": self.n_sims, "draws": self.draws.describe(),
                "screening_start": self.model.screening_start}
        return EnsembleResult(self.model.inputs.unit, self.scenarios, self.model.years,
                              self.model.inputs.denominator, self.model.params,
                              np.ascontiguousarray(C), np.ascontiguousarray(CD),
                              None if full is None else np.ascontiguousarray(full), meta)
