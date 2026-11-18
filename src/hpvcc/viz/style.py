"""
Figure styles of the paper.

Every figure is drawn inside ``with figure_style(name):`` so that styles never leak from
one figure into the next (each figure starts from matplotlib's defaults):

``"paper"``     coverage sweeps and projection plots: the SciencePlots ``science``
                base, the paper's oversized text and light grid, Fira Sans.
``"notebook"``  national-model appendix figures (Figs S3, S4): ``science`` base, 2-pt
                lines, Fira Sans; sizes are set by the figure code itself.
``"plain"``     maps and ranked bar charts: matplotlib defaults with Fira Sans.

The ``science`` settings are copied from SciencePlots 2.x so the package does not need
it installed; LaTeX is never used (all text is Fira Sans with DejaVu Sans maths).
"""
from __future__ import annotations

import contextlib
import os
from functools import lru_cache

import matplotlib

matplotlib.use("Agg")
import matplotlib.font_manager as fm  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from cycler import cycler  # noqa: E402
from matplotlib.ticker import MaxNLocator  # noqa: E402

from ..paths import FONTS_DIR  # noqa: E402

#: SciencePlots "science" style (without its LaTeX settings).
SCIENCE = {
    "axes.linewidth": 0.5,
    "axes.prop_cycle": cycler("color", ["#0C5DA5", "#00B945", "#FF9500", "#FF2C00", "#845B97", "#474747", "#9e9e9e"]),
    "figure.figsize": [3.5, 2.625],
    "grid.linewidth": 0.5,
    "legend.frameon": False,
    "lines.linewidth": 1.0,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.05,
    "xtick.direction": "in", "xtick.major.size": 3.0, "xtick.major.width": 0.5,
    "xtick.minor.size": 1.5, "xtick.minor.visible": True, "xtick.minor.width": 0.5, "xtick.top": True,
    "ytick.direction": "in", "ytick.major.size": 3.0, "ytick.major.width": 0.5,
    "ytick.minor.size": 1.5, "ytick.minor.visible": True, "ytick.minor.width": 0.5, "ytick.right": True,
}

#: Central sizes and colours of the paper figures (oversized: figures are tiled 2–3 per row).
STYLE = dict(
    model_color="#2E7D32", base_color="#D62728",
    model_lw=2.2, grid_alpha=0.22, grid_lw=0.5, grid_nbins=5,
    suptitle_fs=26, title_fs=24, label_fs=24, tick_fs=20, legend_fs=20,
    ci_label_fs=15, baseline_label_fs=18, dpi=130,
)

#: Axis/legend sizes of the national-model figures (Figs 1, S3, S4).
FS_LABEL, FS_TICK, FS_LEGEND, FS_TITLE = 24, 20, 20, 24
GRID_ALPHA, GRID_LW = STYLE["grid_alpha"], STYLE["grid_lw"]

PAPER = {
    "axes.grid": True, "grid.alpha": STYLE["grid_alpha"], "grid.linewidth": STYLE["grid_lw"],
    "axes.titlesize": STYLE["title_fs"], "axes.labelsize": STYLE["label_fs"],
    "xtick.labelsize": STYLE["tick_fs"], "ytick.labelsize": STYLE["tick_fs"],
    "legend.fontsize": STYLE["legend_fs"], "font.size": STYLE["tick_fs"],
    "xtick.minor.visible": False, "ytick.minor.visible": False, "axes.grid.which": "major",
    "legend.frameon": True, "legend.facecolor": "white", "legend.edgecolor": "0.25",
    "legend.framealpha": 1.0, "legend.fancybox": False,
}

FIRA = {"text.usetex": False, "font.family": "Fira Sans", "mathtext.fontset": "dejavusans",
        "axes.unicode_minus": True}

STYLES = {
    "paper": [SCIENCE, PAPER, FIRA],
    "notebook": [SCIENCE, {"lines.linewidth": 2}, FIRA],
    "plain": [{"font.family": "Fira Sans"}],
}


@lru_cache(maxsize=1)
def register_fonts() -> list[str]:
    """Register the bundled Fira Sans fonts with matplotlib (once per process)."""
    added = []
    for f in sorted(FONTS_DIR.glob("*.ttf")):
        fm.fontManager.addfont(str(f))
        added.append(f.name)
    return added


@contextlib.contextmanager
def figure_style(name: str = "paper"):
    """Context manager applying one of :data:`STYLES` on top of matplotlib's defaults."""
    register_fonts()
    with plt.rc_context():
        matplotlib.rcdefaults()
        for layer in STYLES[name]:
            plt.rcParams.update(layer)
        yield


def light_grid(ax, axis: str = "both") -> None:
    """Uniform, sparse major grid (≈5 lines per axis), no minor ticks."""
    ax.grid(True, which="major", axis=axis, alpha=STYLE["grid_alpha"], linewidth=STYLE["grid_lw"])
    ax.minorticks_off()
    ax.yaxis.set_major_locator(MaxNLocator(nbins=STYLE["grid_nbins"]))
    if axis == "both":
        ax.xaxis.set_major_locator(MaxNLocator(nbins=STYLE["grid_nbins"]))


def style_axes(ax, xlabel=None, ylabel=None):
    """Axis labels, tick sizes, ~5 major ticks and a light major grid (national figures)."""
    if xlabel is not None:
        ax.set_xlabel(xlabel, fontsize=FS_LABEL)
    if ylabel is not None:
        ax.set_ylabel(ylabel, fontsize=FS_LABEL)
    ax.tick_params(axis="both", which="major", labelsize=FS_TICK)
    ax.xaxis.set_major_locator(MaxNLocator(5))
    ax.yaxis.set_major_locator(MaxNLocator(5))
    ax.minorticks_off()
    ax.grid(True, which="major", alpha=GRID_ALPHA, linewidth=GRID_LW)


def legend_above(ax, ncol=2, title=None, handles=None, labels=None, fontsize=FS_LEGEND):
    """Boxed legend centred above the axes."""
    box = dict(frameon=True, facecolor="white", edgecolor=(0, 0, 0, 0.25), framealpha=1.0)
    args = () if handles is None else (handles, labels)
    kw = {"title": title, "title_fontsize": FS_TITLE} if title else {}
    return ax.legend(*args, loc="lower center", bbox_to_anchor=(0.5, 1.0), fontsize=fontsize, ncol=ncol, **box, **kw)


# --------------------------------------------------------------------------------------
# Background: figures are saved on white by default; without a background (transparent
# figure and axes, legend boxes kept) on request — the same choice as the model diagrams.
# --------------------------------------------------------------------------------------
_background = {"transparent": os.environ.get("HPVCC_TRANSPARENT", "").lower() in ("1", "true", "yes")}


def set_transparent(on: bool = True) -> None:
    """Save every following figure without (``True``) or with (``False``) a background."""
    _background["transparent"] = bool(on)


def is_transparent() -> bool:
    """Whether figures are currently saved without a background."""
    return _background["transparent"]


@contextlib.contextmanager
def transparent_background(on: bool = True):
    """``with transparent_background(): ...`` — figures saved inside have no background."""
    before = _background["transparent"]
    _background["transparent"] = bool(on)
    try:
        yield
    finally:
        _background["transparent"] = before


def save(fig, path, transparent: bool | None = None, **kwargs):
    """Save a figure (format from the extension of ``path``) and close it.

    ``transparent=None`` follows :func:`set_transparent` / ``HPVCC_TRANSPARENT=1``.
    """
    from pathlib import Path
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    kwargs.setdefault("bbox_inches", "tight")
    if transparent if transparent is not None else _background["transparent"]:
        kwargs["transparent"] = True
    fig.savefig(path, **kwargs)
    plt.close(fig)
    return path
