"""
Per-unit model inputs: projected population, rate denominator, recruitment inflow and
initial conditions.

Inputs are built once from the population and prevalence sources in ``data/raw`` (see
:mod:`hpvcc.data`) and stored as small JSON files in ``data/processed/units/<unit>.json``;
simulations only read those files. They contain no cancer registry data (those are in ``data/registry/``).
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .constants import COMPARTMENTS
from .paths import DATA_PROCESSED
from .units import slug

UNITS_DIR = DATA_PROCESSED / "units"


class RecruitmentInflow:
    """Λ(t): women entering the 18–29 susceptible band per year (thousands).

    Two interpolation rules are used:

    * ``"linear"`` (state models) — linear interpolation between annual values;
    * ``"step"`` (national model) — the value of the calendar year ``floor(t)``,
      clipped to the table's range.
    """

    def __init__(self, years: np.ndarray, values: np.ndarray, mode: str):
        if mode not in ("linear", "step"):
            raise ValueError(f"unknown inflow mode {mode!r}")
        self.years = np.asarray(years, dtype=float)
        self.values = np.asarray(values, dtype=float)
        self.mode = mode
        if mode == "step":
            self._table = {int(y): float(v) for y, v in zip(self.years, self.values)}
            self._lo, self._hi = int(self.years.min()), int(self.years.max())

    def __call__(self, t: float) -> float:
        if self.mode == "linear":
            return float(np.interp(t, self.years, self.values))
        return self._table[int(np.clip(np.floor(t), self._lo, self._hi))]

    def to_dict(self) -> dict:
        return {"years": self.years.tolist(), "values": self.values.tolist(), "mode": self.mode}

    @classmethod
    def from_dict(cls, d: dict) -> "RecruitmentInflow":
        return cls(np.array(d["years"]), np.array(d["values"]), d["mode"])


@dataclass
class UnitInputs:
    """Everything a simulation of one unit needs besides the parameters.

    Attributes
    ----------
    unit : str
        Unit name (``"India"`` or a state/UT).
    years : ndarray
        Simulation years, 2002–2070 (one ODE restart per year).
    population : dict[str, ndarray]
        Projected women by model age band (``"S18"``, ``"S30"``, ``"S61"``), thousands.
    denominator : ndarray
        Women in the rate denominator for each year (thousands). State models use the
        three modelled age bands; the national model uses all-age women.
    inflow : RecruitmentInflow
        Recruitment into the 18–29 band, Λ(t) (before scaling by ``p``).
    y0 : ndarray
        Initial state vector in 2002 (order of :data:`hpvcc.constants.COMPARTMENTS`).
    provenance : dict
        Where each input came from (file names and derivation notes).
    """

    unit: str
    years: np.ndarray
    population: dict
    denominator: np.ndarray
    inflow: RecruitmentInflow
    y0: np.ndarray
    provenance: dict = field(default_factory=dict)

    # ------------------------------------------------------------------ helpers
    def with_y0(self, **values) -> "UnitInputs":
        """Copy with some initial values changed, e.g. ``inputs.with_y0(I=500.0)``."""
        from dataclasses import replace
        y0 = self.y0.copy()
        for name, v in values.items():
            y0[COMPARTMENTS.index(name)] = float(v)
        return replace(self, y0=y0)

    def year_index(self, year: int) -> int:
        return int(np.where(self.years == year)[0][0])

    @property
    def n_years(self) -> int:
        return len(self.years)

    # ------------------------------------------------------------------ persistence
    def to_dict(self) -> dict:
        return {
            "unit": self.unit,
            "years": self.years.tolist(),
            "population": {k: np.asarray(v, dtype=float).tolist() for k, v in self.population.items()},
            "denominator": self.denominator.tolist(),
            "inflow": self.inflow.to_dict(),
            "y0": dict(zip(COMPARTMENTS, self.y0.tolist())),
            "provenance": self.provenance,
        }

    @classmethod
    def from_dict(cls, d: dict) -> "UnitInputs":
        return cls(
            unit=d["unit"],
            years=np.asarray(d["years"], dtype=float),
            population={k: np.asarray(v, dtype=float) for k, v in d["population"].items()},
            denominator=np.asarray(d["denominator"], dtype=float),
            inflow=RecruitmentInflow.from_dict(d["inflow"]),
            y0=np.array([d["y0"][c] for c in COMPARTMENTS], dtype=float),
            provenance=d.get("provenance", {}),
        )

    def save(self, directory: Path | None = None) -> Path:
        directory = Path(directory) if directory else UNITS_DIR
        directory.mkdir(parents=True, exist_ok=True)
        path = directory / f"{slug(self.unit)}.json"
        path.write_text(json.dumps(self.to_dict(), indent=1))
        return path

    @classmethod
    def load(cls, unit: str, directory: Path | None = None) -> "UnitInputs":
        directory = Path(directory) if directory else UNITS_DIR
        path = directory / f"{slug(unit)}.json"
        if not path.exists():
            raise FileNotFoundError(f"{path} not found — run `hpvcc prepare-inputs` first")
        return cls.from_dict(json.loads(path.read_text()))
