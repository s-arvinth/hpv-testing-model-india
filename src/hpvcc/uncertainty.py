"""
Parameter uncertainty of a no-screening projection, from an ensemble of parameter sets.

The national model ships with 500 alternative parameter sets
(``params/india_uncertainty_sets.csv``). They are the "behavioural" sets of a GLUE-style
Monte Carlo analysis done during calibration: 10,000 sets were sampled uniformly within
±25 % of the fitted values of β, β_R, γ, r_I, ε, μ_C, ω and r_C (with r_H ~ U(r_I, 1) and
θ ~ U(0, γ)), scored against the calibration data, and the best 5 % kept, best first.
The scoring step belongs to calibration and is not part of this repository; the accepted
sets are, so the uncertainty band can be regenerated exactly.

>>> ens = ParameterSetEnsemble.from_csv(HPVCervicalCancerModel.for_unit("India"))
>>> res = ens.run()
>>> res.cir_low, res.cir_high          # 95 % band of the no-screening CIR, 2003–2070
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from .metrics import death_counts, incidence_counts
from .model import HPVCervicalCancerModel
from .parameters import ModelParameters
from .paths import PARAMS_DIR

#: Parameters that vary between sets (the others stay at the unit's fitted values).
SAMPLED = ("beta", "beta_reinfection", "gamma", "r_I", "epsilon", "mu_C", "omega", "r_C")

#: Accepted parameter sets of the national model.
INDIA_UNCERTAINTY_SETS = PARAMS_DIR / "india_uncertainty_sets.csv"


@dataclass
class UncertaintyResult:
    """Mean, 2.5/97.5 percentiles and best-set value of CIR and CMR for every year."""

    years: np.ndarray            #: rate years (2003–2070)
    cir_mean: np.ndarray
    cir_low: np.ndarray
    cir_high: np.ndarray
    cmr_mean: np.ndarray
    cmr_low: np.ndarray
    cmr_high: np.ndarray
    cir_best: np.ndarray         #: trajectory of the best-ranked set
    cmr_best: np.ndarray
    sets: np.ndarray             #: (n_sets, len(SAMPLED)) parameter values used
    meta: dict = field(default_factory=dict)

    _ARRAYS = ("years", "cir_mean", "cir_low", "cir_high", "cmr_mean", "cmr_low", "cmr_high",
               "cir_best", "cmr_best", "sets")

    def to_frame(self) -> pd.DataFrame:
        return pd.DataFrame({k: getattr(self, k) for k in self._ARRAYS if k != "sets"})

    def save(self, path: Path) -> Path:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(path, meta=np.array(json.dumps(self.meta)),
                            **{k: getattr(self, k) for k in self._ARRAYS})
        return path

    @classmethod
    def load(cls, path: Path) -> "UncertaintyResult":
        z = np.load(path, allow_pickle=False)
        return cls(**{k: z[k] for k in cls._ARRAYS}, meta=json.loads(str(z["meta"])))


class ParameterSetEnsemble:
    """Simulate the no-screening model once per parameter set and summarise the spread.

    Parameters
    ----------
    model : HPVCervicalCancerModel
        Supplies the inputs and the fitted values of every parameter not in ``sets``.
    sets : DataFrame
        One row per parameter set; columns are parameter names (a subset of
        :class:`~hpvcc.parameters.ModelParameters` fields). Row order = rank (best first).
    """

    def __init__(self, model: HPVCervicalCancerModel, sets: pd.DataFrame, meta: dict | None = None):
        self.model = model
        self.columns = [c for c in sets.columns if c in ModelParameters.__dataclass_fields__]
        self.sets = sets[self.columns].to_numpy(dtype=float)
        self.meta = dict(meta or {})

    @classmethod
    def from_csv(cls, model: HPVCervicalCancerModel, path: Path | None = None, **meta) -> "ParameterSetEnsemble":
        path = Path(path) if path else INDIA_UNCERTAINTY_SETS
        sets = pd.read_csv(path, float_precision="round_trip")
        if "rank" in sets.columns:
            sets = sets.sort_values("rank")
        return cls(model, sets, dict(meta, source=path.name))

    def params_of(self, row: np.ndarray) -> ModelParameters:
        return self.model.params.replace(**dict(zip(self.columns, (float(v) for v in row))))

    def _rates(self, params: ModelParameters):
        # One odeint call over the whole grid.
        traj = self.model.simulate_single_call(params)
        N = np.maximum(self.model.inputs.denominator[1:], 1e-12)
        return incidence_counts(traj[:, 6], params) / N * 1e5, death_counts(traj[:, 9]) / N * 1e5

    def run(self) -> UncertaintyResult:
        cirs, cmrs = [], []
        for row in self.sets:
            try:
                c, m = self._rates(self.params_of(row))
            except Exception:       # a set whose integration fails is skipped
                continue
            cirs.append(c)
            cmrs.append(m)
        cirs, cmrs = np.array(cirs), np.array(cmrs)
        cir_best, cmr_best = self._rates(self.params_of(self.sets[0]))

        def stats(a):
            return np.mean(a, axis=0), np.percentile(a, 2.5, axis=0), np.percentile(a, 97.5, axis=0)

        (cm, cl, ch), (mm, ml, mh) = stats(cirs), stats(cmrs)
        meta = dict(self.meta, n_sets=int(len(self.sets)), n_ok=int(len(cirs)), sampled=self.columns,
                    params=self.model.params.as_dict())
        return UncertaintyResult(self.model.years[1:], cm, cl, ch, mm, ml, mh, cir_best, cmr_best,
                                 self.sets, meta)
