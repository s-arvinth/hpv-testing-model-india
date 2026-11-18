"""
Plots of your own model runs, drawn exactly like the paper's figures.

The bar figures are the paper's own figure functions (:mod:`hpvcc.viz.figures`); the
time-series plots follow the paper's time-series panel (Fig S3, left: SciencePlots base,
Fira Sans, 24-pt labels and titles, 20-pt ticks, 2.5-pt lines with a marker every five
years, 95 % bands at 20 % opacity, boxed legend above the panels).

* :func:`plot_projection` — CIR and CMR over time for every scenario of a
  :class:`~hpvcc.projection.Projection`.
* :func:`plot_coverage_bars` — CIR and CMR by coverage at snapshot years (as Fig 3 / S11–S16).
* :func:`plot_fit` — fitted model vs observed rates, 2012–2019 (as Fig 1 / Fig 2).
* :func:`plot_compartments` — every compartment of a deterministic trajectory.

Every function saves to ``path`` (PDF, PNG or SVG by extension) and returns the path, or
returns the open figure when ``path`` is None. Pass ``transparent=True`` (or call
:func:`hpvcc.viz.style.set_transparent`) to save without a background.
"""
from __future__ import annotations

from pathlib import Path
from typing import Sequence

import matplotlib.pyplot as plt
import numpy as np

from ..constants import COMPARTMENT_LABELS, COMPARTMENT_TEX, COMPARTMENTS, WHO_THRESHOLD
from ..units import INDIA
from .colors import ELIMINATED, INDIA_COLOR, state_color
from .figures import SWEEP_YEARS, TEST_LABELS, coverage_bars, coverage_colors, fit_bars, test_colors
from .style import STYLE, figure_style, save

TITLES = {"CIR": "Crude Incidence Rate (CIR)", "CMR": "Crude Mortality Rate (CMR)"}
BASELINE_LABEL = "No screening (baseline)"
THRESHOLD_LABEL = "CIR = 4 per 100,000 (elimination proxy)"

# Sizes of the paper's time-series panel (Fig S3, left).
_PANEL = (7, 6)
_LW, _MS, _MARK_EVERY, _BAND_ALPHA = 2.5, 6, 5, 0.2
_TITLE_FS, _LABEL_FS, _TICK_FS, _LEGEND_FS = 24, 24, 20, 18


def _out(fig, path, transparent, **kwargs):
    if path is None:
        return fig
    return save(fig, path, transparent=transparent, dpi=kwargs.pop("dpi", 130), **kwargs)


def _legend_above_figure(fig, handles, labels, ncol):
    fig.legend(handles, labels, loc="lower left", bbox_to_anchor=(0.0, 1.0, 1.0, 0.001), mode="expand",
               fontsize=_LEGEND_FS, ncol=ncol, frameon=True, facecolor="white", edgecolor=(0, 0, 0, 0.25),
               framealpha=1.0)


def _style_panel(ax, title=None, xlabel="Year", ylabel=None):
    if title:
        ax.set_title(title, fontsize=_TITLE_FS)
    if xlabel:
        ax.set_xlabel(xlabel, fontsize=_LABEL_FS)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=_LABEL_FS)
    ax.tick_params(axis="both", which="major", labelsize=_TICK_FS)
    ax.grid(True)


def _scenario_label(cov: float, H_pos: float, by_test: bool, by_cov: bool) -> str:
    if cov == 0:
        return BASELINE_LABEL
    test = TEST_LABELS.get(round(H_pos, 2), f"$H^+$ = {H_pos:.2f}")
    parts = []
    if by_cov:
        parts.append(fr"$\alpha$={cov:.2f}")
    if by_test:
        parts.append(test)
    return ", ".join(parts)


def plot_projection(proj, path: str | Path | None = None, metrics: Sequence[str] = ("CIR", "CMR"),
                    start: int = 2012, end: int = 2070, threshold: bool = True,
                    transparent: bool | None = None):
    """CIR and CMR over time, one line (with its 95 % band) per scenario.

    Coverage levels take the unit's colours from the paper (its region's hue family for a
    state/UT, viridis for India); several tests at one coverage take the colours of Fig S3;
    no screening is the red dashed baseline of the coverage figures.

    Parameters
    ----------
    proj : Projection
    path : path, optional
        Save the figure there and close it; otherwise return the open figure.
    metrics : sequence
        ``("CIR",)``, ``("CMR",)`` or both (side by side).
    start, end : int
        Years shown.
    threshold : bool
        Draw the elimination proxy (CIR = 4 per 100,000) on the CIR panel.
    """
    t = proj.table[(proj.table.Year >= start) & (proj.table.Year <= end)]
    groups = list(t.groupby(["Coverage", "Test"], sort=False))
    screened = [(c, s) for (c, s), _ in groups if c > 0]
    covs = sorted({c for c, _ in screened})
    tests = list(dict.fromkeys(s for _, s in screened))
    by_test = len(tests) > 1
    by_cov = len(covs) > 1 or not by_test
    if by_test and not by_cov:
        palette = dict(zip(tests, test_colors(len(tests))))
        color_of = lambda c, s: palette[s]                                  # noqa: E731
    elif by_test:
        palette = dict(zip(screened, test_colors(len(screened))))
        color_of = lambda c, s: palette[(c, s)]                             # noqa: E731
    else:
        levels = coverage_colors(proj.unit, len(covs) + 1)[1:]
        palette = dict(zip(covs, levels))
        color_of = lambda c, s: palette[c]                                  # noqa: E731

    w, h = _PANEL
    with figure_style("notebook"):
        fig, axes = plt.subplots(1, len(metrics), figsize=(w * len(metrics), h), constrained_layout=True,
                                 squeeze=False)
        for ax, metric in zip(axes[0], metrics):
            for (cov, test), g in groups:
                lab = _scenario_label(cov, float(g.H_pos.iloc[0]), by_test, by_cov)
                if cov == 0:
                    ax.plot(g.Year, g[metric], color=STYLE["base_color"], ls="--", lw=_LW, label=lab, zorder=4)
                    continue
                col = color_of(cov, test)
                ax.plot(g.Year, g[metric], color=col, ls="-", lw=_LW, marker="o", markersize=_MS,
                        markevery=_MARK_EVERY, label=lab)
                if (g[f"{metric}_hi95"] > g[f"{metric}_lo95"]).any():
                    ax.fill_between(g.Year, g[f"{metric}_lo95"], g[f"{metric}_hi95"], color=col, alpha=_BAND_ALPHA)
            if metric == "CIR" and threshold:
                ax.axhline(WHO_THRESHOLD, color=ELIMINATED, ls=":", lw=2.0, label=THRESHOLD_LABEL)
            _style_panel(ax, TITLES.get(metric), "Year", f"{metric} / 100,000")
            ax.set_xlim(start, end)
            ax.set_ylim(bottom=0)
        handles, labels = [], []
        for ax in axes[0]:
            for hnd, lab in zip(*ax.get_legend_handles_labels()):
                if lab not in labels:
                    handles.append(hnd)
                    labels.append(lab)
        _legend_above_figure(fig, handles, labels, ncol=min(len(labels), 4 if len(metrics) > 1 else 2))
        return _out(fig, path, transparent)


def plot_coverage_bars(proj, path: str | Path | None = None, years: Sequence[int] = SWEEP_YEARS,
                       transparent: bool | None = None):
    """CIR and CMR by coverage at a few years, exactly as Fig 3 (India) / Figs S11–S16 (states).

    ``proj`` must hold one screening test; no screening is added if it was not projected.
    With the paper's coverages (0, 0.25, 0.5, 0.7, 0.9, 1) and the default ensemble, a
    state/UT figure is identical to its panel in the paper.
    """
    snap, covs = proj.coverage_snapshots(years)
    return coverage_bars(proj.unit, snap, path, years, covs, transparent=transparent)


def plot_fit(unit: str, path: str | Path | None = None, model=None, transparent: bool | None = None):
    """Fitted (no-screening) model vs observed CIR and CMR, 2012–2019 — Fig 1 for India,
    Fig 2 / Figs S5–S10 for a state/UT. ``model`` may be a modified model."""
    from ..observed import model_fit, observed_label
    f = model_fit(unit, model)
    unit = f.Unit.iloc[0]
    return fit_bars(unit, f.Year, f.Model_CIR.to_numpy(), f.Observed_CIR.to_numpy(), f.Model_CMR.to_numpy(),
                    f.Observed_CMR.to_numpy(), path, observed_label=observed_label(unit), transparent=transparent)


def plot_compartments(model, traj: np.ndarray | None = None, path: str | Path | None = None,
                      start: int | None = None, transparent: bool | None = None):
    """Every compartment of a deterministic trajectory (thousands of women), one panel each,
    in the unit's colour."""
    traj = model.simulate() if traj is None else traj
    years = model.years
    keep = np.ones(len(years), bool) if start is None else years >= start
    color = INDIA_COLOR if model.unit == INDIA else state_color(model.unit)
    w, h = _PANEL
    with figure_style("notebook"):
        fig, axes = plt.subplots(2, 5, figsize=(5 * w, 2 * h), constrained_layout=True)
        for k, (ax, name) in enumerate(zip(axes.flat, COMPARTMENTS)):
            ax.plot(years[keep], traj[keep, k], color=color, ls="-", lw=_LW, marker="o", markersize=_MS,
                    markevery=_MARK_EVERY)
            _style_panel(ax, f"{COMPARTMENT_LABELS[name]} (${COMPARTMENT_TEX[name]}$)",
                         "Year" if k >= 5 else None, "Population (10$^3$)" if k % 5 == 0 else None)
        return _out(fig, path, transparent)


__all__ = ["plot_projection", "plot_coverage_bars", "plot_fit", "plot_compartments"]
