#!/usr/bin/env python3
"""Plot reference science spectra and transmittance by profile and geometry."""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.cm as cm
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import netCDF4 as nc
import numpy as np

from validationlib import load_cases

ROOT = Path(__file__).resolve().parents[1]
REFS = ROOT / "references"
PLOTS = ROOT / "plots"
GEOMETRIES = ("zenith", "nadir", "limb")
COLORMAP = LinearSegmentedColormap.from_list(
    "geometry", ["darkblue", "blue", "cyan", "green", "yellow", "orange", "red"]
)


def variable_wavenumber(name: str) -> float:
    return float(name.split("_", 1)[1])


def load_case(case_name: str, geometry: str, prefix: str, profile: str):
    path = REFS / case_name / "rad.nc"
    if not path.is_file():
        raise FileNotFoundError(f"missing {profile} reference: {path}")
    with nc.Dataset(path) as dataset:
        names = sorted((name for name in dataset.variables if name.startswith(prefix)), key=variable_wavenumber)
        if not names:
            raise ValueError(f"no {prefix} variables in {path}")
        values = np.stack([
            np.ma.filled(dataset.variables[name][0, :], np.nan).astype(float) for name in names
        ], axis=1)
        if geometry == "limb":
            geometry_values = np.ma.filled(dataset.variables["tp_z"][0, :], np.nan).astype(float)
            geometry_label = "Tangent altitude (km)"
        else:
            geometry_values = np.ma.filled(dataset.variables["vp_lat"][0, :], np.nan).astype(float)
            geometry_label = "View-point latitude (deg)"
    return np.array([variable_wavenumber(name) for name in names]), values, geometry_values, geometry_label


def collect_spectrum(geometry: str, prefix: str, profile: str):
    cases = [case for case in load_cases() if case.profile == profile and case.geometry == geometry]
    cases.sort(key=lambda case: (case.nu_min, case.nu_max))
    spectra: dict[float, np.ndarray] = {}
    geometry_values = None
    geometry_label = ""
    for case in cases:
        nus, values, case_geometry, geometry_label = load_case(
            case.case_name, geometry, prefix, profile
        )
        if geometry_values is None:
            geometry_values = case_geometry
        elif not np.allclose(geometry_values, case_geometry, equal_nan=True):
            raise ValueError(f"geometry differs between {profile} {geometry} cases")
        for index, nu in enumerate(nus):
            if nu in spectra and not np.allclose(spectra[nu], values[:, index], equal_nan=True):
                raise ValueError(f"overlapping chunks disagree at {nu:.4f} cm^-1")
            spectra.setdefault(nu, values[:, index])
    if not spectra or geometry_values is None:
        raise ValueError(f"no {profile} cases found for geometry: {geometry}")
    nus = np.array(sorted(spectra))
    values = np.stack([spectra[nu] for nu in nus], axis=1)
    return cases, nus, values, geometry_values, geometry_label


def plot_dataset(geometry: str, kind: str, profile: str) -> None:
    science_prefix = "rad_" if geometry == "limb" else "bt_"
    prefix = "tau_" if kind == "tau" else science_prefix
    cases, nus, values, geometry_values, geometry_label = collect_spectrum(
        geometry, prefix, profile
    )
    finite_geometry = geometry_values[np.isfinite(geometry_values)]
    geo_min, geo_max = float(np.min(finite_geometry)), float(np.max(finite_geometry))
    if geo_min == geo_max:
        geo_max = geo_min + 1.0
    norm = plt.Normalize(geo_min, geo_max)
    fig, ax = plt.subplots(figsize=(14, 6), facecolor="white", constrained_layout=True)
    for ray_index in range(values.shape[0]):
        color = COLORMAP(norm(geometry_values[ray_index]))
        ax.plot(nus, values[ray_index], color=color, linewidth=0.6, alpha=0.6)
    ax.set_xlabel("Wavenumber (cm$^{-1}$)")
    if kind == "tau":
        ax.set_ylabel("Transmittance")
        ax.set_ylim(-0.02, 1.02)
        quantity = "Transmittance"
    elif geometry == "limb":
        ax.set_ylabel("Radiance (W m$^{-2}$ sr$^{-1}$ (cm$^{-1}$)$^{-1}$)")
        ax.set_yscale("log")
        quantity = "Radiance"
    else:
        ax.set_ylabel("Brightness temperature (K)")
        quantity = "Brightness Temperature"
    ax.set_title(
        f"{profile.capitalize()} Reference {quantity}: {geometry.capitalize()}\n"
        f"{len(cases)} cases, {nus[0]:.0f}–{nus[-1]:.0f} cm$^{{-1}}$, {values.shape[0]} rays"
    )
    ax.grid(alpha=0.3, linewidth=0.5)
    ax.set_xlim(nus[0] - 5, nus[-1] + 5)
    colorbar = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=COLORMAP), ax=ax, shrink=0.8)
    colorbar.set_label(geometry_label)
    PLOTS.mkdir(parents=True, exist_ok=True)
    suffix = "_tau" if kind == "tau" else ""
    output = PLOTS / f"reference_{profile}_{geometry}{suffix}.png"
    fig.savefig(output, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {output}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--profile",
        choices=("full", "smoke"),
        default="full",
        help="validation profile to plot (default: full)",
    )
    args = parser.parse_args()

    for selected_geometry in GEOMETRIES:
        plot_dataset(selected_geometry, "science", args.profile)
        plot_dataset(selected_geometry, "tau", args.profile)


if __name__ == "__main__":
    main()
