"""
Fitted parameter sets.

The final parameters are stored as plain CSV in ``params/``:

* ``params/states_final.csv`` — one row per state/UT model (36 rows);
* ``params/india_final.csv`` — the national model (one row).

Both files hold full double-precision values so that simulations reproduce the
paper's numbers exactly. To run the model with your own values, either change them in
Python (``params.replace(beta=0.4)``) or write a CSV in the same format and pass it to
:func:`load_parameters` (``file=...``) or to the command line (``hpvcc run --params-file``).
Columns that a file leaves out take the defaults μ = 0.01388 and H⁺ = 0.95.

The clearance rate r_H and the progression rate θ of screen-detected infection are not
stored: they are sampled in every simulation (r_H ~ U(r_I, 1), θ ~ U(0, γ), redrawn every
year; see :mod:`hpvcc.simulation`).
"""
from __future__ import annotations

import dataclasses
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from .constants import H_POS_HPV_DNA, MU_NATURAL
from .paths import PARAMS_DIR
from .units import INDIA, STATES

@dataclass(frozen=True)
class ModelParameters:
    """Parameters of the 10-compartment model for one unit.

    Rates are per year. Populations are in thousands of women. The two rates of the
    screened pathway, r_H and θ, are sampled and therefore not part of a parameter set
    (:data:`SAMPLED_RATES`).

    Attributes
    ----------
    p : float
        Scaling of the recruitment inflow Λ(t) into the 18–29 susceptible band.
    K1, K2 : float
        Ageing rates 18–29 → 30–65 and 30–65 → 66–85.
    beta : float
        HPV transmission rate (identical in the three susceptible age bands).
    beta_reinfection : float
        Transmission rate for re-susceptible women (S_R).
    gamma : float
        Progression rate from undetected infection (I) to cancer (C).
    r_I : float
        Clearance rate of undetected infection.
    r_C : float
        Rate of leaving the cancer compartment to the survivor compartment.
    mu_C, n_PL : float
        Cancer mortality, ``mu_C * C**n_PL`` (power law).
    omega : float
        Relapse rate from survivors (C_S) back to cancer.
    epsilon : float
        Rate of loss of post-clearance immunity (R_H → S_R).
    mu : float
        Natural mortality.
    H_pos : float
        Screening-test sensitivity (probability an infected screened woman tests positive).
    """

    p: float
    K1: float
    K2: float
    beta: float
    beta_reinfection: float
    gamma: float
    r_I: float
    r_C: float
    mu_C: float
    n_PL: float
    omega: float
    epsilon: float
    mu: float = MU_NATURAL
    H_pos: float = H_POS_HPV_DNA

    def replace(self, **changes) -> "ModelParameters":
        """Return a copy with some fields changed (e.g. ``params.replace(H_pos=0.70)``)."""
        return dataclasses.replace(self, **changes)

    def as_dict(self) -> dict[str, float]:
        return dataclasses.asdict(self)


#: Symbol, meaning and unit of every parameter (used by the docs and ``hpvcc params``).
PARAMETER_INFO: dict[str, tuple[str, str, str]] = {
    "p": ("p", "scaling of the recruitment inflow Λ(t) into S18", "–"),
    "K1": ("K1", "ageing 18–29 → 30–65", "per year"),
    "K2": ("K2", "ageing 30–65 → 66–85", "per year"),
    "beta": ("β", "HPV transmission rate (S18, S30, S61)", "per year"),
    "beta_reinfection": ("β_R", "transmission rate for re-susceptible women (S_R)", "per year"),
    "gamma": ("γ", "progression of undetected infection to cancer (I → C)", "per year"),
    "r_I": ("r_I", "clearance of undetected infection (I → R_H)", "per year"),
    "r_C": ("r_C", "cancer → survivor (C → C_S)", "per year"),
    "mu_C": ("μ_C", "cancer mortality coefficient, deaths = μ_C·C^n_PL", "–"),
    "n_PL": ("n_PL", "exponent of the cancer-mortality power law", "–"),
    "omega": ("ω", "relapse of survivors (C_S → C)", "per year"),
    "epsilon": ("ε", "loss of post-clearance immunity (R_H → S_R)", "per year"),
    "mu": ("μ", "natural mortality, all compartments", "per year"),
    "H_pos": ("H⁺", "screening-test sensitivity", "probability"),
}

#: The sampled rates of the screened pathway: symbol, meaning, distribution.
SAMPLED_RATES: dict[str, tuple[str, str, str]] = {
    "r_H": ("r_H", "clearance / treatment of screen-detected infection (H_T → R_H)",
            "U(r_I, 1), redrawn every year"),
    "theta": ("θ", "progression of screen-detected infection to cancer (H_T → C)",
              "U(0, γ), redrawn every year"),
}

# --------------------------------------------------------------------------------------
# Loaders
# --------------------------------------------------------------------------------------

STATE_PARAMS_CSV = PARAMS_DIR / "states_final.csv"
INDIA_PARAMS_CSV = PARAMS_DIR / "india_final.csv"

#: CSV column → field name, where they differ.
_COLUMN_ALIASES = {"beta_const": "beta", "beta_Reinfec": "beta_reinfection"}


def params_from_row(row: pd.Series) -> ModelParameters:
    """Parameters from one CSV row (column names as in ``params/*.csv`` or field names).

    Missing ``mu`` and ``H_pos`` take the defaults (μ = 0.01388, H⁺ = 0.95). Columns
    ``theta`` and ``r_H`` are ignored: these rates are sampled.
    """
    d = {_COLUMN_ALIASES.get(k, k): v for k, v in row.items()}

    def get(name, default=None):
        v = d.get(name, default)
        if v is None or (isinstance(v, float) and pd.isna(v)):
            if default is None:
                raise KeyError(f"parameter {name!r} missing")
            return float(default)
        return float(v)

    return ModelParameters(
        p=get("p"), K1=get("K1"), K2=get("K2"), beta=get("beta"),
        beta_reinfection=get("beta_reinfection"), gamma=get("gamma"), r_I=get("r_I"),
        r_C=get("r_C"), mu_C=get("mu_C"), n_PL=get("n_PL"), omega=get("omega"),
        epsilon=get("epsilon"), mu=get("mu", MU_NATURAL), H_pos=get("H_pos", H_POS_HPV_DNA),
    )


def _read(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, float_precision="round_trip")


def load_parameters(unit: str, params_dir: Path | None = None, file: Path | None = None) -> ModelParameters:
    """Fitted parameters of ``unit`` (``"India"`` or a state/UT name).

    ``file`` reads a CSV of your own instead: either one row (used for any unit) or one
    row per unit with a ``State`` (or ``Unit``) column.
    """
    from .units import resolve_unit
    unit = resolve_unit(unit)
    if file is not None:
        df = _read(Path(file))
        key = next((c for c in ("State", "Unit") if c in df.columns), None)
        if key is None:
            return params_from_row(df.iloc[0])
        df = df.set_index(key)
        if unit not in df.index:
            raise KeyError(f"no row for {unit!r} in {file}")
        return params_from_row(df.loc[unit])
    params_dir = Path(params_dir) if params_dir else PARAMS_DIR
    if unit == INDIA:
        return params_from_row(_read(params_dir / INDIA_PARAMS_CSV.name).iloc[0])
    df = _read(params_dir / STATE_PARAMS_CSV.name).set_index("State")
    if unit not in df.index:
        raise KeyError(f"No parameters for {unit!r} in {STATE_PARAMS_CSV.name}")
    return params_from_row(df.loc[unit])


def parameters_frame(units: list[str] | None = None) -> pd.DataFrame:
    """Fitted parameters of several units as one table (rows = units, columns = fields)."""
    units = units or [INDIA, *STATES]
    return pd.DataFrame([load_parameters(u).as_dict() for u in units], index=units)


def load_all_parameters(params_dir: Path | None = None) -> dict[str, ModelParameters]:
    """Parameters of India and all 36 state/UT models."""
    return {u: load_parameters(u, params_dir) for u in [INDIA, *STATES]}
