"""
State choropleth maps of projected CIR (Fig 4).

One fixed colour scale per metric across all maps (0 → the highest state value over the
reference grid of coverages 0/25/50/75/90/100 % and years 2030–2070), so maps are directly
comparable. States whose mean CIR is below 4 per 100,000 are hatched.
Requires ``geopandas`` (``pip install hpvcc[maps]``).
"""
from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Patch

from ...constants import WHO_THRESHOLD
from ..experiments import Experiment
from ...paths import GEO_DIR
from ...viz.colors import EMBER, darken
from ...viz.style import figure_style, save

ELIMINATION_LABEL = "CIR < 4 per 100,000 (elimination proxy)"
#: Grid over which the fixed colour-scale maximum is taken.
SCALE_COVERAGES = (0.0, 0.25, 0.50, 0.75, 0.90, 1.00)
SCALE_YEARS = tuple(range(2030, 2071, 5))


def _clip(name: str, lo: float, hi: float) -> LinearSegmentedColormap:
    c = matplotlib.colormaps[name]
    return LinearSegmentedColormap.from_list(f"{name}_c", c(np.linspace(hi, lo, 256)))


#: Colour schemes (the paper uses ``02_magma_noblack``).
CMAPS = {
    "01_YlOrRd": "YlOrRd", "02_magma_noblack": ("magma", 0.24, 0.96), "03_inferno_noblack": ("inferno", 0.22, 0.96),
    "04_ember_custom": EMBER, "05_plasma": "plasma", "06_viridis": "viridis", "07_cividis": "cividis",
    "08_Reds": "Reds", "09_YlGnBu": "YlGnBu", "10_Spectral_r": "Spectral_r", "11_RdYlGn_r": "RdYlGn_r",
    "12_turbo": "turbo",
}


def colormap(scheme: str):
    spec = CMAPS[scheme]
    if isinstance(spec, tuple):
        return _clip(*spec)
    return matplotlib.colormaps[spec] if isinstance(spec, str) else spec


@lru_cache(maxsize=1)
def states_geodataframe():
    import geopandas as gpd
    gdf = gpd.read_file(GEO_DIR / "india_states_dissolved.geojson").to_crs(7755)
    return gdf.set_index("st_nm").loc[sorted(gdf["st_nm"])].reset_index()   # alphabetical drawing order


@lru_cache(maxsize=None)
def state_values(metric: str, year: int, coverage: float) -> dict[str, float]:
    """Mean CIR or CMR of every state in ``year`` at ``coverage``."""
    exp = Experiment.get("states_coverage_sweep")
    out = {}
    for unit in exp.units:
        res = exp.result(unit)
        s = res.snapshot(res.index(coverage=coverage), year)
        out[unit] = s.cir if metric == "CIR" else s.cmr
    return out


@lru_cache(maxsize=None)
def scale_max(metric: str) -> float:
    return float(max(v for Y in SCALE_YEARS for a in SCALE_COVERAGES for v in state_values(metric, Y, a).values()))


def choropleth(metric: str, year: int, coverage: float, out: Path, scheme: str = "02_magma_noblack",
               title: bool = False) -> Path:
    gdf = states_geodataframe()
    states = list(gdf["st_nm"])
    vals = state_values(metric, year, coverage)
    cir = state_values("CIR", year, coverage)
    arr = np.array([vals[s] for s in states])
    elim = np.array([cir[s] < WHO_THRESHOLD for s in states])
    cmap = colormap(scheme)
    norm = Normalize(0, scale_max(metric))
    with figure_style("plain"):
        face = [cmap(norm(v)) for v in arr]
        edge = [darken(c) for c in face]
        fig, ax = plt.subplots(figsize=(9.5, 11))
        ax.set_axis_off()
        gdf.plot(ax=ax, color=face, edgecolor=edge, linewidth=0.9)
        if elim.any():
            sub = gdf[elim]
            ec = [darken(cmap(norm(arr[i]))) for i in np.where(elim)[0]]
            sub.plot(ax=ax, facecolor="none", edgecolor=ec, hatch="///", linewidth=1.6)
        sm = cm.ScalarMappable(cmap=cmap, norm=norm)
        cb = fig.colorbar(sm, ax=ax, orientation="horizontal", fraction=0.045, pad=0.02, shrink=0.62)
        cb.set_label(f"{metric} per 100,000", fontsize=12)
        cb.ax.tick_params(labelsize=10)
        ax.legend(handles=[Patch(facecolor="none", edgecolor="0.25", hatch="///", label=ELIMINATION_LABEL)],
                  loc="lower left", fontsize=10, frameon=True)
        if title:
            name = {"CIR": "Crude Incidence Rate", "CMR": "Crude Mortality Rate"}[metric]
            fig.suptitle(f"{name} (per 100,000) — {year}, {round(coverage * 100)}% screening coverage",
                         fontsize=14, fontweight="bold", y=0.965)
        return save(fig, out)
