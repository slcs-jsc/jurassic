"""Publication-style scientific PNG/SVG plots with a headless Matplotlib backend.

This module renders existing JURASSIC results; it never runs or modifies the C model.
"""
from __future__ import annotations

from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

from .parsers import chart, atmospheric_chart, kernel_heatmap

COLORS = ("#12bec8", "#8656e8", "#ef8c47", "#49aa76", "#eb5690", "#156fae")


def _style(fig, axes, theme):
    if theme not in ("light", "dark"):
        raise ValueError("Invalid plot theme")
    is_dark = theme == "dark"
    background = "#182338" if is_dark else "#ffffff"
    panel = "#202e45" if is_dark else "#ffffff"
    foreground = "#e7effb" if is_dark else "#22314c"
    grid = "#54637b" if is_dark else "#b3c0d2"
    fig.patch.set_facecolor(background)
    for ax in axes:
        ax.set_facecolor(panel)
        ax.tick_params(colors=foreground, labelsize=9)
        ax.xaxis.label.set_color(foreground)
        ax.yaxis.label.set_color(foreground)
        ax.title.set_color(foreground)
        for spine in ax.spines.values():
            spine.set_color(grid)
        ax.grid(alpha=0.21, color=grid)
    return foreground


def _save(fig, output):
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)


def render_chart(source, params, metric, output, theme="light"):
    """Export forward-model brightness temperature/radiance or transmittance."""
    info = chart(source, params["geometry"], params["channels"], metric)
    fig, ax = plt.subplots(figsize=(10, 5.2), dpi=150)
    fg = _style(fig, (ax,), theme)
    for i, trace in enumerate(info["traces"]):
        ax.plot(info["x"], trace["y"], label=trace["name"], linewidth=2.3,
                color=COLORS[i % len(COLORS)])
    ax.set(xlabel=info["xlabel"], ylabel=info["ylabel"])
    ax.set_title(f"JURASSIC · {params['geometry'].title()} · " +
                 ("Transmittance" if metric == "transmittance" else
                  "Brightness Temperature" if params.get("output") == "bbt" else "Radiance"))
    legend = ax.legend(fontsize=9, facecolor=fig.get_facecolor(), framealpha=0.86)
    for label in legend.get_texts(): label.set_color(fg)
    fig.tight_layout()
    _save(fig, output)


def render_atmosphere(source, gases, output, theme="light"):
    """Display temperature, pressure, and gas mixing ratios on separate scales."""
    data = atmospheric_chart(source, gases)
    fig, axes = plt.subplots(1, 3, figsize=(12.6, 5.3), dpi=145, sharey=True)
    fg = _style(fig, axes, theme)
    temp = next(t for t in data["traces"] if t["type"] == "temperature")
    pressure = next(t for t in data["traces"] if t["type"] == "pressure")
    axes[0].plot(temp["x"], temp["y"], color=COLORS[0], linewidth=2.4)
    axes[0].set(xlabel="Temperature (K)", ylabel="Altitude (km)", title="Temperature")
    axes[1].plot(pressure["x"], pressure["y"], color=COLORS[1], linewidth=2.4)
    if all(v > 0 for v in pressure["x"]): axes[1].set_xscale("log")
    axes[1].set(xlabel="Pressure (hPa)", title="Pressure")
    gases_data = [t for t in data["traces"] if t["type"] == "gas"]
    for i, trace in enumerate(gases_data):
        axes[2].plot(trace["x"], trace["y"], label=trace["name"].replace(" mixing ratio", ""),
                     color=COLORS[i % len(COLORS)], linewidth=2)
    if gases_data and all(v > 0 for t in gases_data for v in t["x"]):
        axes[2].set_xscale("log")
    axes[2].set(xlabel="Volume mixing ratio (ppv)", title="Selected gases")
    if gases_data:
        legend = axes[2].legend(loc="best", fontsize=8, facecolor=fig.get_facecolor(), framealpha=0.85)
        for label in legend.get_texts(): label.set_color(fg)
    fig.suptitle("JURASSIC · Atmospheric profile", color=fg, fontsize=13, fontweight="bold")
    fig.tight_layout()
    _save(fig, output)


def render_kernel(source, output, theme="light"):
    """Render a sampled temperature Jacobian heatmap from an existing kernel.tab."""
    data = kernel_heatmap(source, "TEMPERATURE", max_rows=110, max_altitudes=150)
    z = np.array([[np.nan if v is None else float(v) for v in row] for row in data["z"]])
    finite = z[np.isfinite(z)]
    if finite.size == 0:
        raise ValueError("Kernel contains no finite temperature sensitivities")
    magnitude = float(np.max(np.abs(finite))) or 1.0
    cmap = plt.get_cmap("coolwarm").copy()
    cmap.set_bad("#78879a" if theme == "dark" else "#ecf2f9")
    fig, ax = plt.subplots(figsize=(10, 5.2), dpi=150)
    fg = _style(fig, (ax,), theme)
    x0, x1 = float(data["x"][0]), float(data["x"][-1])
    y0, y1 = int(data["y"][0]), int(data["y"][-1])
    if x0 == x1: x0 -= 0.5; x1 += 0.5
    if y0 == y1: y0 -= 0.5; y1 += 0.5
    im = ax.imshow(np.ma.masked_invalid(z), origin="lower", aspect="auto", interpolation="nearest",
                   extent=(x0, x1, y0, y1), cmap=cmap,
                   norm=TwoSlopeNorm(vmin=-magnitude, vcenter=0, vmax=magnitude))
    ax.set(xlabel="Altitude (km)", ylabel="Measurement index",
           title="Temperature sensitivity · first 110 measurement indices")
    cbar = fig.colorbar(im, ax=ax, fraction=0.047, pad=0.025)
    cbar.set_label("Jacobian element", color=fg)
    cbar.ax.tick_params(colors=fg)
    cbar.outline.set_edgecolor(fg)
    fig.tight_layout()
    _save(fig, output)
