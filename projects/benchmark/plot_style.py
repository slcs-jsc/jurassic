"""Shared plot style for the benchmark evaluation scripts.

Publication-style defaults: a clean sans-serif face, recessive hairline grid,
only left/bottom spines, thin 2 px lines with ringed markers, and a
colourblind-checked categorical palette (first three slots validated for all
pairs, first eight for adjacent pairs). Text always wears the ink colours, never
a series colour.

Usage:
    from plot_style import apply_style, INK, PALETTE, save_figure
    apply_style()
"""
from __future__ import annotations

import textwrap
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt

# --- colour tokens (light surface) ---
SURFACE = "#ffffff"
INK = "#0b0b0b"          # primary text
INK_2 = "#52514e"        # secondary text, axis labels
INK_MUTED = "#898781"    # reference lines (ideal scaling, ceilings)
GRID = "#e6e5e0"         # hairline grid
AXIS = "#b9b8b2"         # spines and ticks

# Categorical slots in fixed order (never cycled): blue, orange, aqua, yellow,
# magenta, green, violet, red.
PALETTE = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100",
           "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
BLUE, ORANGE, AQUA = PALETTE[0], PALETTE[1], PALETTE[2]

FONT_STACK = ["Helvetica", "Arial", "Liberation Sans", "Lato", "DejaVu Sans"]

# Figure sizes [inch]: fits a single column (~3.4) or a full text width (~6.5)
FIG_SINGLE = (6.4, 4.0)


def apply_style() -> None:
    """Install the rcParams used by all benchmark plots."""
    mpl.rcParams.update({
        # canvas
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "figure.figsize": FIG_SINGLE,
        "figure.dpi": 110,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        # type
        "font.family": "sans-serif",
        "font.sans-serif": FONT_STACK,
        "mathtext.fontset": "stixsans",
        "font.size": 10,
        "axes.titlesize": 11,
        "axes.titleweight": "semibold",
        "axes.titlelocation": "left",
        "axes.titlepad": 10,
        "axes.labelsize": 10,
        "axes.labelcolor": INK_2,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.fontsize": 9,
        "text.color": INK,
        "xtick.color": AXIS,
        "ytick.color": AXIS,
        "xtick.labelcolor": INK_2,
        "ytick.labelcolor": INK_2,
        # axes: only left and bottom, thin
        "axes.edgecolor": AXIS,
        "axes.linewidth": 0.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "xtick.major.size": 3.5,
        "ytick.major.size": 3.5,
        "xtick.minor.size": 2.0,
        "ytick.minor.size": 2.0,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        # recessive solid hairline grid, behind the data
        "axes.grid": True,
        "axes.grid.axis": "y",
        "axes.axisbelow": True,
        "grid.color": GRID,
        "grid.linewidth": 0.8,
        "grid.linestyle": "-",
        # data marks: 2 px lines, ringed markers
        "lines.linewidth": 1.6,
        "lines.markersize": 6.5,
        "lines.markeredgewidth": 1.2,
        "lines.solid_capstyle": "round",
        "lines.solid_joinstyle": "round",
        "axes.prop_cycle": mpl.cycler(color=PALETTE),
        # legend: no frame, no shadow
        "legend.frameon": False,
        "legend.handlelength": 1.8,
        "legend.borderaxespad": 0.4,
        "legend.labelcolor": INK_2,
    })


def wrap(text: str, width: int = 34) -> str:
    """Wrap a long axis label onto several lines."""
    return "\n".join(textwrap.wrap(text, width=width, break_long_words=False))


def line_marker_kwargs(color: str) -> dict:
    """Series line with a filled marker ringed in the surface colour."""
    return dict(color=color, marker="o", markerfacecolor=color,
                markeredgecolor=SURFACE)


def save_figure(fig, path: Path, formats=("png", "pdf")) -> None:
    """Save `fig` to `path`; also write the other listed formats next to it
    (vector PDF for papers). `path` may end in .png or have no suffix."""
    path = Path(path)
    stem = path.with_suffix("")
    for fmt in formats:
        fig.savefig(f"{stem}.{fmt}")
    plt.close(fig)
