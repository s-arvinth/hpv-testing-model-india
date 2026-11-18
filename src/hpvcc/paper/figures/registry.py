"""
Registry of every figure in the paper and its appendix.

Each figure is a list of panels; each panel is one output file, written under
``results/paper/figures/`` with the same relative path as in the manuscript's figure
folder, so the output can be dropped into the manuscript unchanged.

Figures 1, 2 and S5–S10 compare the fitted models with the observed registry rates
(``data/registry/``); every other figure is drawn from the paper's simulations.

>>> from hpvcc.paper.figures.registry import make_figure
>>> make_figure("fig3")        # runs the India coverage sweep if needed, then draws Fig 3
"""
from __future__ import annotations

from dataclasses import dataclass, field
from functools import partial
from pathlib import Path
from typing import Callable

import pandas as pd

from ...paths import FIGURES_DIR, FIGURES_TRANSPARENT_DIR, TABLES_DIR
from ...units import REGIONS, slug
from .calibration import fit_figure
from .national import model_diagram, monte_carlo_bars, test_type_figure
from .sweeps import coverage_sweep_bars

PAPER_FIGURES = FIGURES_DIR

@dataclass
class Panel:
    path: str                               #: output path relative to the paper figure folder
    draw: Callable[[Path], Path]            #: callable(out_path) -> out_path


@dataclass
class PaperFigure:
    key: str
    number: str
    label: str
    section: str
    caption: str
    experiments: list[str]
    panels: list[Panel] = field(default_factory=list)


def _annual_long() -> pd.DataFrame:
    """Master annual table (built and cached on first use)."""
    from ...paths import SIMULATIONS_DIR
    csv = TABLES_DIR / "csv" / "Annual_Long.csv"
    sims = [p for d in ("states_coverage_sweep", "india_coverage_sweep") for p in (SIMULATIONS_DIR / d).glob("*.npz")]
    stale = csv.exists() and sims and max(p.stat().st_mtime for p in sims) > csv.stat().st_mtime
    if not csv.exists() or stale:
        from ..tables import MasterTables
        MasterTables().save()
    return pd.read_csv(csv)


def _ranked(metric: str, year: int, cov: int, out: Path) -> Path:
    from .ranked import figure_source, ranked_bars
    return ranked_bars(figure_source(_annual_long()), metric, year, cov, out, layout="grid")


def _map(year: int, cov: float, out: Path) -> Path:
    from .maps import choropleth
    return choropleth("CIR", year, cov, out, scheme="02_magma_noblack")


_CAL = "State_Level_Figures/Calibration_Bars_Fira/Rate_Calibration_Bars/{}_rate_calib.pdf"
_SWEEP = "State_Level_Figures/AlphaSweep_YearCombos_Fira_70pct/2040_2050_2060/{}_alpha_sweep_2040_2050_2060.pdf"


def _state_cal(unit):
    return Panel(_CAL.format(slug(unit)), partial(fit_figure, unit))


def _state_sweep(unit):
    return Panel(_SWEEP.format(slug(unit)), partial(coverage_sweep_bars, unit))


def _region_members(region):
    return sorted(REGIONS[region])


def build_registry() -> dict[str, PaperFigure]:
    F = {}

    def add(key, number, label, section, caption, experiments, panels):
        F[key.lower()] = PaperFigure(key.lower(), number, label, section, caption, experiments, panels)

    add("fig1", "Fig 1", "fig:national_bar_calib", "main",
        "National model vs registry CIR and CMR, 2012-2019", [],
        [Panel("India_Level_Figures/Calibration_Bars/CIR_CMR_bar_comparison_calib_window.pdf",
               partial(fit_figure, "India"))])
    add("fig2", "Fig 2", "fig:region_calib_representatives", "main",
        "Calibration of one representative state per region", [],
        [_state_cal(u) for u in ("Madhya Pradesh", "West Bengal", "Delhi", "Arunachal Pradesh",
                                 "Tamil Nadu", "Maharashtra")])
    add("fig3", "Fig 3", "fig:national_alpha_snapshot", "main",
        "National CIR and CMR by screening coverage, 2040/2050/2060", ["india_coverage_sweep"],
        [Panel("India_Level_Figures/India_AlphaSweep_YearCombos_Fira_70pct_viridis/India_alpha_sweep_2040_2050_2060.pdf",
               partial(coverage_sweep_bars, "India"))])
    add("fig4", "Fig 4", "fig:cir_map_grid", "main",
        "State CIR maps at 0/25/70/90 % coverage, 2040/2050/2060", ["states_coverage_sweep"],
        [Panel(f"India_Maps/Publication_Maps/02_magma_noblack/untitled/CIR/coverage_{c:03d}pct/CIR_{y}.pdf",
               partial(_map, y, c / 100)) for c in (0, 25, 70, 90) for y in (2040, 2050, 2060)])
    for key, number, label, metric in (("fig5", "Fig 5", "fig:cir_rankhist_grid", "CIR"),
                                       ("fig6", "Fig 6", "fig:cmr_rankhist_grid", "CMR")):
        add(key, number, label, "main", f"Ranked {metric} of all units at 0 % and 70 % coverage, 2040/2050",
            ["states_coverage_sweep", "india_coverage_sweep"],
            [Panel(f"Region_Comparisons/Portrait/Ranked_Histograms_Grid/03_elimination_status/untitled/{metric}/"
                   f"coverage_{c:03d}pct/{metric}_{y}.pdf", partial(_ranked, metric, y, c))
             for y in (2040, 2050) for c in (0, 70)])
    add("figS1", "Fig S1", "fig:finaltesting", "appendix", "Model diagram with HPV testing", [],
        [Panel("Model_Figure/HPV_Testing_Model_for_Paper.pdf", partial(model_diagram, "full"))])
    add("figS2", "Fig S2", "fig:baselinemodel", "appendix", "Baseline model diagram", [],
        [Panel("Model_Figure/HPV_Testing_Baseline_Model_for_Paper.pdf", partial(model_diagram, "baseline"))])
    add("figS3", "Fig S3", "fig:national_hpos", "appendix",
        "National model: four screening tests at 75 % coverage", ["india_test_type"],
        [Panel("India_Level_Figures/Modality_Comparison_India/GoldStandard_Hpos_combined_detection_snapshot.pdf",
               test_type_figure)])
    add("figS4", "Fig S4", "fig:mc_cir", "appendix",
        "Monte Carlo uncertainty of the national no-screening projection", ["india_monte_carlo"],
        [Panel(f"India_Level_Figures/Monte_Carlo/MonteCarlo_{m}_barchart_sparse.pdf", partial(monte_carlo_bars, m))
         for m in ("CIR", "CMR")])
    regions = [("Central", "central"), ("East", "east"), ("North", "north"), ("Northeast", "northeast"),
               ("South", "south"), ("West", "west")]
    for i, (region, tag) in enumerate(regions):
        add(f"figS{5 + i}", f"Fig S{5 + i}", f"fig:app_region_calib_{tag}", "appendix",
            f"Calibration, {region} region", [], [_state_cal(u) for u in _region_members(region)])
    suffix = {"North": "_a", "Northeast": "_a", "South": "_a"}
    for i, (region, tag) in enumerate(regions):
        add(f"figS{11 + i}", f"Fig S{11 + i}", f"fig:app_region_alpha_{tag}{suffix.get(region, '')}", "appendix",
            f"Coverage sweep, {region} region", ["states_coverage_sweep"],
            [_state_sweep(u) for u in _region_members(region)])
    return F


REGISTRY = build_registry()


def make_figure(key: str, out_dir: Path | None = None, verbose: bool = True,
                transparent: bool = False) -> list[Path]:
    """Draw every panel of one paper figure (``key`` like ``"fig3"``, ``"figS11"`` or ``"Fig S11"``).

    ``transparent=True`` saves the panels without a background (and the white-line model
    diagrams) under ``results/paper/figures_transparent/`` unless ``out_dir`` is given.
    """
    from ...viz.style import transparent_background
    norm = key.lower().replace(" ", "")
    if norm not in REGISTRY:
        raise KeyError(f"unknown figure {key!r}; choose from {', '.join(REGISTRY)}")
    fig = REGISTRY[norm]
    out_dir = Path(out_dir) if out_dir else (FIGURES_TRANSPARENT_DIR if transparent else PAPER_FIGURES)
    paths = []
    with transparent_background(transparent):
        for panel in fig.panels:
            paths.append(panel.draw(out_dir / panel.path))
            if verbose:
                print(f"  {fig.number}: {panel.path}")
    return paths


def make_all(out_dir: Path | None = None, verbose: bool = True, transparent: bool = False) -> dict[str, list[Path]]:
    """Draw every figure of the paper and its appendix."""
    return {k: make_figure(k, out_dir, verbose, transparent) for k in REGISTRY}
