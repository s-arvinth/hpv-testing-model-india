"""
Ranked horizontal bars of CIR or CMR across the 36 states/UTs and India (Figs 5–6).

Bars are sorted by value, with the 95 % CI as error bars and India (dedicated national
model) highlighted. Colour = elimination status: green when the unit's mean CIR is below
4 per 100,000 in that year and coverage; otherwise (CIR figure) an ember scale from 4 to
the panel maximum, or (CMR figure) grey. Values are taken from the master table rounded
to three decimals, as in the paper.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize
from matplotlib.patches import Patch

from ...constants import WHO_THRESHOLD
from ...units import INDIA
from ...viz.colors import ELIMINATED, EMBER, INDIA_BAND, INDIA_COLOR, NOT_ELIMINATED, darken
from ...viz.style import figure_style, save

LABEL_BELOW = "CIR < 4 per 100,000 (elimination proxy)"
LABEL_ABOVE = "CIR ≥ 4 per 100,000"

LAYOUTS = {
    "single": dict(fig=(9, 12.4), ytick=10, ytick_in=11, xtick=11, xlab=12.5, val=9, val_in=10.5, leg=9.5, elw=1.0, cap=2.5),
    "grid": dict(fig=(6.4, 7.9), ytick=10.5, ytick_in=11.5, xtick=10.5, xlab=11.5, val=9, val_in=10.5, leg=9.5, elw=0.8,
                 cap=1.8, leg_out=True),
}


def figure_source(annual_long: pd.DataFrame) -> pd.DataFrame:
    """CIR/CMR mean and CI half-width per unit, year and coverage, rounded as in the paper."""
    cols = {"Unit": "State", "CIR_mean": "CIR", "CMR_mean": "CMR", "CIR_hw95": "CIR_95CI_hw", "CMR_hw95": "CMR_95CI_hw"}
    df = annual_long[list(cols) + ["Year", "Coverage_pct"]].rename(columns=cols)
    for c in ("CIR", "CMR", "CIR_95CI_hw", "CMR_95CI_hw"):
        df[c] = df[c].round(3)
    return df


def ranked_bars(source: pd.DataFrame, metric: str, year: int, coverage_pct: int, out: Path,
                layout: str = "grid", title: bool = False) -> Path:
    L = LAYOUTS[layout]
    errkw = dict(ecolor="0.25", elinewidth=L["elw"], capsize=L["cap"])
    sub = source[(source.Year == year) & (source.Coverage_pct == coverage_pct)].set_index("State")
    states = sorted(s for s in sub.index if s != INDIA)
    units = states + [INDIA]
    vals = np.array([float(sub.loc[u, metric]) for u in units])
    cis = np.array([float(sub.loc[u, f"{metric}_95CI_hw"]) for u in units])
    elim = np.array([float(sub.loc[u, "CIR"]) < WHO_THRESHOLD for u in units])
    order = np.argsort(vals)
    units = [units[i] for i in order]
    vals, cis, elim = vals[order], cis[order], elim[order]
    is_india = np.array([u == INDIA for u in units])
    n_top, i_india = len(units) - 1, int(np.where(is_india)[0][0])

    if metric == "CIR":
        cols = [ELIMINATED if v < WHO_THRESHOLD else EMBER(Normalize(WHO_THRESHOLD, max(vals))(v)) for v in vals]
        legend = [Patch(fc=ELIMINATED, label=LABEL_BELOW), Patch(fc="#f16913", label=LABEL_ABOVE)]
    else:
        cols = [ELIMINATED if e else NOT_ELIMINATED for e in elim]
        legend = [Patch(fc=ELIMINATED, label=LABEL_BELOW), Patch(fc=NOT_ELIMINATED, label=LABEL_ABOVE)]

    with figure_style("plain"):
        xm = (vals + cis).max() * 1.14
        fig, ax = plt.subplots(figsize=L["fig"])
        ax.set_xlabel(f"{metric} per 100,000  (mean, 95% CI)", fontsize=L["xlab"])
        ax.set_xlim(0, xm)
        ax.set_ylim(-0.8, n_top + 2.2)
        ax.set_yticks(range(len(units)))
        ax.set_yticklabels(units, fontsize=L["ytick"])
        for t, ind in zip(ax.get_yticklabels(), is_india):
            if ind:
                t.set_fontweight("bold"); t.set_color(INDIA_COLOR); t.set_fontsize(L["ytick_in"])
        ax.tick_params(axis="x", labelsize=L["xtick"])
        ax.grid(axis="x", alpha=0.25, lw=0.6)
        ax.set_axisbelow(True)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.axhspan(i_india - 0.5, i_india + 0.5, color=INDIA_BAND, zorder=0)
        bars = ax.barh(range(len(units)), vals, color=cols, edgecolor=[darken(c) for c in cols], linewidth=0.9,
                       xerr=cis, error_kw=errkw)
        bars[i_india].set_edgecolor(INDIA_COLOR); bars[i_india].set_linewidth(2.6); bars[i_india].set_zorder(4)
        for i, (v, e, ind) in enumerate(zip(vals, cis, is_india)):
            ax.text(v + e + xm * 0.008, i, f"{v:.1f}", va="center", ha="left",
                    fontsize=L["val_in"] if ind else L["val"], fontweight="bold" if ind else "normal",
                    color=INDIA_COLOR if ind else "0.15")
        if L.get("leg_out"):
            ax.legend(handles=legend, loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2, fontsize=L["leg"],
                      frameon=True, handlelength=1.4, columnspacing=1.2, borderaxespad=0.3)
        else:
            ax.legend(handles=legend, loc="lower right", fontsize=L["leg"], frameon=True)
        if title:
            name = {"CIR": "Crude Incidence Rate", "CMR": "Crude Mortality Rate"}[metric]
            fig.suptitle(f"{name} (per 100,000) — {year}, {coverage_pct}% screening coverage",
                         fontsize=14, fontweight="bold", y=0.995)
        fig.tight_layout(rect=[0, 0, 1, 0.985])
        return save(fig, out)
