"""
The figure types of the paper, as functions of data.

Every figure of the manuscript is drawn by one of these functions, and so is every figure
the package draws for you (:meth:`hpvcc.projection.Projection.plot_bars`,
:func:`hpvcc.viz.plots.plot_fit`, ...): same layout, fonts, sizes, colours, legends and
ticks. The paper subpackage (:mod:`hpvcc.paper.figures`) only collects the paper's
simulations and passes them in.

=========================  =====================================  ======================
Function                   Figure type                            In the paper
=========================  =====================================  ======================
:func:`coverage_bars`      CIR/CMR by coverage at snapshot years  Fig 3, Figs S11–S16
:func:`fit_bars`           fitted model vs observed, 2012–2019    Fig 1 (India), Fig 2,
                                                                  Figs S5–S10 (states)
:func:`test_bars`          four screening tests                   Fig S3
:func:`uncertainty_bars`   parameter-uncertainty band             Fig S4
=========================  =====================================  ======================

Each function saves to ``out`` (PDF, PNG or SVG by extension) and returns the path; with
``out=None`` it returns the open figure instead.
"""
from __future__ import annotations

from pathlib import Path
from typing import Mapping, Sequence

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import AutoMinorLocator, MaxNLocator

from ..units import INDIA
from .colors import complement, region_cmap, region_of, state_color
from .style import (GRID_ALPHA, GRID_LW, STYLE, figure_style, legend_above, light_grid, save,
                    style_axes)

#: Coverages and snapshot years of the coverage figures (0.70 = WHO target).
SWEEP_COVERAGES = (0.0, 0.25, 0.50, 0.70, 0.90, 1.00)
SWEEP_YEARS = (2040, 2050, 2060)

#: Legend labels of the four screening tests (by sensitivity).
TEST_LABELS = {0.70: "Pap Cytology (0.70)", 0.75: "VIA (0.75)", 0.85: "LBC (0.85)", 0.95: "HPV DNA Test (0.95)"}

MC_COLORS = {"CIR": "#B5651D", "CMR": "#37474F"}
MC_TITLES = {"CIR": "Crude Incidence Rate (CIR)", "CMR": "Crude Mortality Rate (CMR)"}


def _finish(fig, out, transparent=None, **kwargs):
    """Save inside the style context (savefig settings are part of the style)."""
    if out is None:
        return fig
    return save(fig, out, transparent=transparent, **kwargs)


# ============================================================================ colours
def viridis_sweep_cmap() -> LinearSegmentedColormap:
    """Viridis remapped so that sampling it at 0.95 → 0.30 spans viridis(0) → viridis(1)."""
    xs = np.linspace(0, 1, 256)
    return LinearSegmentedColormap.from_list("viridis_combo", plt.cm.viridis(np.clip((0.95 - xs) / 0.65, 0, 1)))


def coverage_colors(unit: str, n: int) -> np.ndarray:
    """Colours of ``n`` coverage levels (no screening first = darkest): the unit's region hue
    family for a state/UT, remapped viridis for India."""
    if unit == INDIA:
        cmap = viridis_sweep_cmap()
    else:
        reg = region_of(unit)
        cmap = region_cmap(reg) if reg else plt.cm.Greens
    return cmap(np.linspace(0.95, 0.30, n))


def test_colors(n: int) -> np.ndarray:
    """Colours of ``n`` screening tests, lowest sensitivity first."""
    return plt.cm.cool(np.linspace(0, 1, n))


# ============================================================================ coverage
def coverage_bars(unit: str, snap: Mapping[float, Mapping[int, Sequence[float]]], out: Path | None = None,
                  years: Sequence[int] = SWEEP_YEARS, coverages: Sequence[float] = SWEEP_COVERAGES,
                  transparent: bool | None = None):
    """CIR (top) and CMR (bottom) by screening coverage at a few snapshot years (Fig 3, S11–S16).

    ``snap[coverage][year] = (CIR, CIR half-width, CMR, CMR half-width)``. The first coverage
    is the reference drawn as the red dashed "No screening" segment, so it should be 0.
    One bar per coverage (mean ± 95 % half-width, value printed on top).
    """
    colors = coverage_colors(unit, len(coverages))
    red = STYLE["base_color"]
    n_a = len(coverages)
    bw = 0.10
    sp = n_a * bw + 0.45
    gc = [i * sp for i in range(len(years))]
    off = [(j - n_a / 2 + 0.5) * bw for j in range(n_a)]
    with figure_style("paper"):
        fig, axs = plt.subplots(2, 1, figsize=(14, 10), constrained_layout=True)
        for ai, (a, col) in enumerate(zip(coverages, colors)):
            xp = [gc[gi] + off[ai] for gi in range(len(years))]
            for ax, mi in ((axs[0], 0), (axs[1], 2)):
                m = [snap[a][Y][mi] for Y in years]
                e = [snap[a][Y][mi + 1] for Y in years]
                ax.bar(xp, m, bw, color=col, alpha=0.9, yerr=e, capsize=3,
                       error_kw={"lw": 1, "ecolor": "k"}, label=fr"$\alpha$={a:.2f}")
                for x_, m_, e_ in zip(xp, m, e):
                    ax.text(x_, m_ + e_, f"{m_:.1f}", ha="center", va="bottom", fontsize=STYLE["ci_label_fs"],
                            rotation=90, color="0.15")
        for ax, mi, ylabel in ((axs[0], 0, "CIR / 100,000"), (axs[1], 2, "CMR / 100,000")):
            for gi, Y in enumerate(years):
                base = snap[coverages[0]][Y][mi]
                x0, x1 = gc[gi] + off[0] - bw * 0.9, gc[gi] + off[-1] + bw * 0.9
                ax.hlines(base, x0, x1, color=red, ls="--", lw=1.8, zorder=7)
                ax.annotate(f"{base:.1f}", (gc[gi], base), xytext=(0, 4), textcoords="offset points",
                            fontsize=STYLE["baseline_label_fs"], color=red, ha="center", va="bottom",
                            fontweight="bold")
            ax.plot([], [], color=red, ls="--", lw=1.8, label="No screening (baseline)")
            ax.set_xticks(gc)
            ax.set_xticklabels([str(Y) for Y in years])
            ax.set_xlabel("Year")
            ax.set_ylabel(ylabel)
            ax.set_ylim(0, ax.get_ylim()[1] * 1.22)
            ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=7, columnspacing=1.0, handletextpad=0.4)
            light_grid(ax, "y")
        return _finish(fig, out, transparent, dpi=STYLE["dpi"])


# ============================================================================ model fit
def fit_bars(unit: str, years, cir_model, cir_observed, cmr_model, cmr_observed, out: Path | None = None,
             observed_label: str = "Observed", transparent: bool | None = None):
    """Fitted model vs observed CIR (left) and CMR (right) per year, 2012–2019.

    India is drawn as Fig 1 (blue model bars, red observed bars); a state/UT as Figs 2 and
    S5–S10 (the state's own colour vs the complementary hue).
    """
    years = np.asarray(years).astype(int)
    if unit == INDIA:
        return _national_fit_bars(years, cir_model, cir_observed, cmr_model, cmr_observed, out, observed_label,
                                  transparent)
    with figure_style("paper"):
        fig, axs = plt.subplots(1, 2, figsize=(14, 5.2), constrained_layout=True)
        x = np.arange(len(years))
        bw = 0.38
        for ax, mv, ov, ylabel in ((axs[0], cir_model, cir_observed, "CIR / 100,000"),
                                   (axs[1], cmr_model, cmr_observed, "CMR / 100,000")):
            ax.bar(x - bw / 2, mv, bw, color=state_color(unit), alpha=0.9, label="Model")
            ax.bar(x + bw / 2, ov, bw, color=complement(unit), alpha=0.9, label=observed_label)
            ax.set_xticks(x)
            ax.set_xticklabels(years, rotation=45, ha="right")
            ax.set_ylabel(ylabel)
            ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2, borderaxespad=0.2)
            light_grid(ax, "y")
        return _finish(fig, out, transparent, dpi=STYLE["dpi"])


def _national_fit_bars(years, cir_m, cir_o, cmr_m, cmr_o, out, observed_label, transparent=None):
    with figure_style("notebook"):
        fig, axs = plt.subplots(1, 2, figsize=(14, 5.5), constrained_layout=True)
        bw = 0.35
        x = np.arange(len(years))
        for ax, mv, ov, mc, oc, ml, ylabel in (
                (axs[0], cir_m, cir_o, "royalblue", "firebrick", "Model CIR", "CIR / 100,000"),
                (axs[1], cmr_m, cmr_o, "cornflowerblue", "darkred", "Model CMR", "CMR / 100,000")):
            ax.bar(x - bw / 2, mv, bw, color=mc, alpha=0.85, label=ml)
            ax.bar(x + bw / 2, ov, bw, color=oc, alpha=0.85, label=observed_label)
            style_axes(ax, ylabel=ylabel)
            ax.set_xticks(x)
            ax.set_xticklabels(years, rotation=45, ha="right", fontsize=STYLE["tick_fs"])
            legend_above(ax, ncol=2)
        return _finish(fig, out, transparent, dpi=130)


# ============================================================================ screening tests
def test_bars(years_axis, detected, snaps, sensitivities, out: Path | None = None,
              years: Sequence[int] = (2030, 2050, 2070), transparent: bool | None = None):
    """Screen-detected infections over time (left) and CIR / CMR at a few years (right) for
    several screening tests at one coverage (Fig S3).

    Parameters
    ----------
    years_axis : array
        Years of the trajectories.
    detected : list of TrajectoryStats
        ``H_T`` (mean, 2.5 %, 97.5 %) per test.
    snaps : list of dict
        ``snaps[i][year]`` = :class:`~hpvcc.metrics.RateSnapshot` of test ``i``.
    sensitivities : list of float
    """
    colors = test_colors(len(sensitivities))
    bar_colors = colors
    labels = [TEST_LABELS.get(round(h, 2), f"$H^+$ = {h:.2f}") for h in sensitivities]

    left_w, left_h, right_w, right_h, height = 7, 6, 13, 10, 6
    right_ratio = left_w * (right_w / right_h) / (left_w / left_h)
    with figure_style("notebook"):
        fig = plt.figure(figsize=(left_w + right_ratio, height), constrained_layout=True)
        fig.get_layout_engine().set(h_pad=0.01, w_pad=0.04)
        gs = fig.add_gridspec(nrows=2, ncols=2, width_ratios=[left_w, right_ratio], height_ratios=[1, 1])
        ax_l, ax_cir, ax_cmr = fig.add_subplot(gs[:, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, 1])

        handles = []
        for st, col in zip(detected, colors):
            line, = ax_l.plot(years_axis, st.mean, color=col, linestyle="-", linewidth=2.5, marker="o",
                              markersize=6, markevery=5)
            ax_l.fill_between(years_axis, st.low, st.high, color=col, alpha=0.2)
            handles.append(line)
        ax_l.set_title(r"HPV Detection ($H_T$)", fontsize=24)
        ax_l.set_xlabel("Year", fontsize=24)
        ax_l.set_ylabel("Population (10$^3$)", fontsize=24)
        ax_l.tick_params(axis="both", which="major", labelsize=20)
        ax_l.grid(True)

        n_h, bw = len(sensitivities), 0.12
        group = n_h * bw + 0.45
        centres = [i * group for i in range(len(years))]
        offsets = [(j - n_h / 2 + 0.5) * bw for j in range(n_h)]
        for si, col in enumerate(bar_colors):
            s_y = [snaps[si][Y] for Y in years]
            xpos = [centres[g] + offsets[si] for g in range(len(years))]
            for ax, m, e in ((ax_cir, [s.cir for s in s_y], [s.cir_halfwidth for s in s_y]),
                             (ax_cmr, [s.cmr for s in s_y], [s.cmr_halfwidth for s in s_y])):
                ax.bar(xpos, m, bw, color=col, alpha=0.85, yerr=e, capsize=3, error_kw={"lw": 1.0, "ecolor": "black"})
                for x, v, err in zip(xpos, m, e):
                    ax.text(x, v + err, f"{v:.2f}", ha="center", va="bottom", fontsize=15, rotation=90)
        for ax, ylabel in ((ax_cir, "CIR / 100,000"), (ax_cmr, "CMR / 100,000")):
            style_axes(ax, xlabel="Year", ylabel=ylabel)
            ax.set_xlabel("Year", fontsize=20)
            ax.set_ylabel(ylabel, fontsize=20)
            ax.tick_params(axis="both", which="major", labelsize=18)
            ax.set_xticks(centres)
            ax.set_xticklabels([str(Y) for Y in years], fontsize=18)
            ax.set_ylim(top=ax.get_ylim()[1] * 1.30)
        fig.legend(handles, labels, loc="lower left", bbox_to_anchor=(0.0, 1.0, 1.0, 0.001), mode="expand",
                   fontsize=18, ncol=4, frameon=True, facecolor="white", edgecolor=(0, 0, 0, 0.25), framealpha=1.0)
        return _finish(fig, out, transparent, dpi=130)


# ============================================================================ parameter uncertainty
def uncertainty_bars(result, metric: str, out: Path | None = None, step: int = 10, width: float = 6.0,
                     n_label: int | None = None, transparent: bool | None = None):
    """Best-fit value with the parameter-uncertainty 95 % band every ``step`` years (Fig S4).

    ``result`` is an :class:`~hpvcc.uncertainty.UncertaintyResult`; ``n_label`` is the
    number of sampled parameter sets quoted in the legend (default: the sets used).
    """
    low, high, best = {"CIR": (result.cir_low, result.cir_high, result.cir_best),
                       "CMR": (result.cmr_low, result.cmr_high, result.cmr_best)}[metric]
    yrs = result.years
    bar_years = np.arange(int(np.ceil(yrs[0])), int(np.floor(yrs[-1])) + 1, step)
    idx = [int(np.argmin(np.abs(yrs - Y))) for Y in bar_years]
    x = yrs[idx]
    ticks = np.arange(int(np.ceil(x[0] / 10) * 10), int(np.floor(x[-1] / 10) * 10) + 1, 10)
    color = MC_COLORS[metric]
    with figure_style("notebook"):
        fig, ax = plt.subplots(figsize=(7.5, 7.2), constrained_layout=True)
        b, lo, hi = best[idx], low[idx], high[idx]
        bars = ax.bar(x, b, width=width, color=color, alpha=0.85, yerr=np.vstack([b - lo, hi - b]), capsize=7,
                      error_kw={"lw": 1.8, "ecolor": color, "alpha": 1.0, "capthick": 2.0}, zorder=2)
        style_axes(ax, xlabel="Year", ylabel=f"{metric} / 100,000")
        ax.set_xticks(list(ticks))
        ax.set_xticklabels([str(int(t)) for t in ticks], rotation=0, ha="center")
        ax.yaxis.set_major_locator(MaxNLocator(8))
        ax.minorticks_on()
        ax.yaxis.set_minor_locator(AutoMinorLocator(4))
        ax.grid(True, which="major", axis="y", alpha=GRID_ALPHA, linewidth=GRID_LW)
        ax.grid(True, which="minor", axis="y", alpha=GRID_ALPHA * 0.5, linewidth=GRID_LW * 0.7)
        ax.set_axisbelow(True)
        patch = mpatches.Patch(facecolor=color, alpha=0.85, edgecolor=color, linewidth=1.0)
        n = n_label if n_label is not None else len(result.sets)
        legend_above(ax, ncol=1, title=MC_TITLES[metric], handles=[patch, bars.errorbar],
                     labels=["Best-fit value", f"95% Confidence Band (Monte Carlo, $n={n}$)"])
        return _finish(fig, out, transparent, dpi=130)


__all__ = ["coverage_bars", "fit_bars", "test_bars", "uncertainty_bars", "coverage_colors", "test_colors",
           "viridis_sweep_cmap", "SWEEP_COVERAGES", "SWEEP_YEARS", "TEST_LABELS"]
