"""
One colour per state/UT, used in every figure.

Each region owns a hue family; within a region, states are shaded by their mean model
(no-screening) CIR over 2012–2019 (lowest = lightest, highest = darkest).
"""
from __future__ import annotations

import colorsys
from functools import lru_cache

import numpy as np
from matplotlib.colors import LinearSegmentedColormap, to_hex, to_rgb

from ..units import REGIONS, REGION_OF

#: Region → (light end, dark end) of its hue family.
REGION_ENDS = {
    "South": ("#cfe3f7", "#08306b"),      # blues
    "North": ("#fdd0c0", "#7f0000"),      # reds
    "Northeast": ("#d3efcf", "#00441b"),  # greens
    "West": ("#e3d5f0", "#3f007d"),       # purples
    "East": ("#fddcb0", "#7f2704"),       # oranges
    "Central": ("#c6ece8", "#00504a"),    # teal
}
_LO, _HI = 0.42, 0.92          # shade range used within a hue family

INDIA_COLOR = "#7a0177"
INDIA_BAND = "#f3e0ef"
ELIMINATED = "#1a9850"
NOT_ELIMINATED = "#4d4d4d"
EMBER = LinearSegmentedColormap.from_list(
    "ember", ["#fff7ec", "#fdd49e", "#fdae61", "#f16913", "#d7301f", "#7f0000"])


def region_cmap(region: str) -> LinearSegmentedColormap:
    lo, hi = REGION_ENDS[region]
    return LinearSegmentedColormap.from_list(region, [lo, hi])


@lru_cache(maxsize=1)
def baseline_cir_2012_2019() -> dict[str, float]:
    """Mean model CIR 2012–2019 of every state (from the no-screening trajectories)."""
    from ..metrics import rates
    from ..model import HPVCervicalCancerModel
    out = {}
    for members in REGIONS.values():
        for unit in members:
            m = HPVCervicalCancerModel.for_unit(unit)
            cir, _ = rates(m.simulate(), m.params, m.inputs.denominator)
            yrs = list(m.years[1:].astype(int))
            out[unit] = float(np.nanmean([cir[yrs.index(y)] for y in range(2012, 2020)]))
    return out


@lru_cache(maxsize=1)
def state_colors() -> dict[str, str]:
    """Fixed colour of every state/UT."""
    cir = baseline_cir_2012_2019()
    colors = {}
    for region, members in REGIONS.items():
        ordered = sorted(members, key=lambda s: (np.isnan(cir[s]), cir[s]))
        cm = region_cmap(region)
        pos = np.linspace(_LO, _HI, len(ordered)) if len(ordered) > 1 else np.array([0.75])
        colors.update({s: to_hex(cm(p)) for s, p in zip(ordered, pos)})
    return colors


def state_color(unit: str) -> str:
    return state_colors().get(unit, "#2E7D32")


def complement(unit: str) -> str:
    """Colour of the observed series beside ``unit``'s model series (opposite hue)."""
    r, g, b = to_rgb(state_color(unit))
    h, s, v = colorsys.rgb_to_hsv(r, g, b)
    return to_hex(colorsys.hsv_to_rgb((h + 0.5) % 1.0, 0.62, 0.78))


def darken(color, k: float = 0.55):
    """Same hue, lower lightness — used for bar and polygon edges."""
    r, g, b = to_rgb(color)
    h, l, s = colorsys.rgb_to_hls(r, g, b)
    return colorsys.hls_to_rgb(h, max(0, l * k), min(1, s * 1.05))


def region_of(unit: str) -> str | None:
    return REGION_OF.get(unit)
