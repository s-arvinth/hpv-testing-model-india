"""
The named simulations behind the paper (``paper/experiments.yaml``): run once, cache, reuse.

>>> exp = Experiment.get("states_coverage_sweep")
>>> res = exp.result("Kerala")            # loads the cached ensemble, or runs it first
>>> res.snapshot(res.index(coverage=0.70), 2050).cir
"""
from __future__ import annotations

import time
from functools import lru_cache
from pathlib import Path

import yaml

from ..model import HPVCervicalCancerModel
from ..paths import EXPERIMENTS_FILE, ROOT, SIMULATIONS_DIR
from ..simulation import (EnsembleResult, StochasticEnsemble, coverage_sweep, default_workers,
                          make_draw_scheme, sensitivity_sweep)
from ..uncertainty import ParameterSetEnsemble, UncertaintyResult
from ..units import STATES, slug


@lru_cache(maxsize=1)
def load_config(path: str | None = None) -> dict:
    with open(path or EXPERIMENTS_FILE) as fh:
        return yaml.safe_load(fh)["experiments"]


class Experiment:
    """One named experiment from ``paper/experiments.yaml``.

    Two kinds exist: ``ensemble`` (stochastic screening scenarios, the default) and
    ``parameter_sets`` (the national parameter-uncertainty band of Fig S4).
    """

    def __init__(self, name: str, spec: dict):
        self.name = name
        self.spec = spec

    @classmethod
    def get(cls, name: str) -> "Experiment":
        cfg = load_config()
        if name not in cfg:
            raise KeyError(f"unknown experiment {name!r}; choose from {sorted(cfg)}")
        return cls(name, cfg[name])

    @classmethod
    def all(cls) -> list["Experiment"]:
        return [cls(n, s) for n, s in load_config().items()]

    # ------------------------------------------------------------------ properties
    @property
    def description(self) -> str:
        return " ".join(self.spec.get("description", "").split())

    @property
    def units(self) -> list[str]:
        u = self.spec["units"]
        return list(STATES) if u == "states" else list(u)

    @property
    def is_parameter_sets(self) -> bool:
        return self.spec.get("kind") == "parameter_sets"

    def scenarios(self):
        values, fixed = self.spec["values"], self.spec.get("fixed", {})
        if self.spec["sweep"] == "coverage":
            return coverage_sweep(values, fixed["H_pos"])
        return sensitivity_sweep(values, fixed["coverage"])

    def path(self, unit: str) -> Path:
        return SIMULATIONS_DIR / self.name / f"{slug(unit)}.npz"

    # ------------------------------------------------------------------ run / load
    def run(self, unit: str, n_sims: int | None = None, workers: int | None = None,
            save: bool = True, verbose: bool = True):
        """Run the experiment for one unit (``n_sims`` overrides the configured number)."""
        model = HPVCervicalCancerModel.for_unit(unit)
        t0 = time.time()
        if self.is_parameter_sets:
            res = ParameterSetEnsemble.from_csv(model, ROOT / self.spec["sets"],
                                                n_samples_drawn=self.spec["n_samples_drawn"]).run()
        else:
            res = StochasticEnsemble(model, self.scenarios(), n_sims or self.spec["n_sims"],
                                     make_draw_scheme(self.spec["draws"], self.scenarios()),
                                     keep_full=self.spec.get("keep_full", False),
                                     workers=workers or default_workers()).run()
            res.meta["experiment"] = self.name
        if verbose:
            print(f"  [{self.name}] {unit}: {time.time() - t0:.0f} s")
        if save:
            res.save(self.path(unit))
        return res

    def load(self, unit: str):
        loader = UncertaintyResult.load if self.is_parameter_sets else EnsembleResult.load
        return loader(self.path(unit))

    def _is_current(self, unit: str, res) -> bool:
        """Was ``res`` made with the current parameters and experiment size?"""
        from ..parameters import load_parameters
        if self.is_parameter_sets:
            import numpy as np
            import pandas as pd
            sets = pd.read_csv(ROOT / self.spec["sets"], float_precision="round_trip")
            return (res.meta.get("params") == load_parameters(unit).as_dict()
                    and np.array_equal(res.sets, sets[res.meta["sampled"]].to_numpy(dtype=float)))
        return (res.params == load_parameters(unit) and res.n_sims == self.spec["n_sims"]
                and len(res.scenarios) == len(self.spec["values"]))

    def result(self, unit: str, **run_kwargs):
        """Cached result, (re-)running the experiment first if it is missing or out of date
        (parameters edited in ``params/`` or experiment size changed in the config)."""
        if self.path(unit).exists():
            res = self.load(unit)
            if self._is_current(unit, res):
                return res
        return self.run(unit, **run_kwargs)

    def run_all(self, **kwargs):
        for unit in self.units:
            self.run(unit, **kwargs)
