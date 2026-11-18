"""
Raw data → model inputs.

``data/raw/`` holds the source files (see ``data/raw/SOURCES.md``); the builders here turn
them into one :class:`hpvcc.inputs.UnitInputs` per unit, written to ``data/processed/units``:

* :mod:`hpvcc.data.states`   — the 36 state/UT models
* :mod:`hpvcc.data.national` — the national (India) model
"""
from __future__ import annotations

from ..units import INDIA, STATES


def build_all(verbose: bool = True):
    """Build and save the inputs of every unit. Returns ``{unit: path}``."""
    from .national import build_national_inputs
    from .sources import RawData
    from .states import build_state_inputs

    raw = RawData()
    paths = {INDIA: build_national_inputs(raw).save()}
    for unit in STATES:
        paths[unit] = build_state_inputs(unit, raw).save()
        if verbose:
            print(f"  inputs: {unit}")
    return paths
