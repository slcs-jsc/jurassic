#!/usr/bin/env python3
"""Plot per-channel science and transmittance differences for one validation case."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np


def wavenumber(name: str) -> float:
    return float(name.split("_", 1)[1])


def error_series(reference, candidate, prefix: str):
    ref_names = sorted((name for name in reference.variables if name.startswith(prefix)), key=wavenumber)
    test_names = sorted((name for name in candidate.variables if name.startswith(prefix)), key=wavenumber)
    if ref_names != test_names:
        raise ValueError(
            f"reference and candidate {prefix} variable sets differ "
            f"(reference={len(ref_names)}, candidate={len(test_names)}); "
            "the run is likely stale relative to the current references"
        )
    nus, absolute, relative = [], [], []
    for name in ref_names:
        ref = np.ma.filled(reference.variables[name][:], np.nan).astype(float)
        test = np.ma.filled(candidate.variables[name][:], np.nan).astype(float)
        if ref.shape != test.shape:
            raise ValueError(f"shape mismatch for {name}: {ref.shape} != {test.shape}")
        ref_finite = np.isfinite(ref)
        test_finite = np.isfinite(test)
        if not np.array_equal(ref_finite, test_finite):
            raise ValueError(f"finite/NaN mask mismatch for {name}")
        finite = ref_finite & test_finite
        diff = np.abs(test[finite] - ref[finite])
        ref_values = ref[finite]
        nonzero = ref_values != 0
        nus.append(wavenumber(name))
        absolute.append(float(np.max(diff)) if diff.size else np.nan)
        relative.append(
            float(np.max(100.0 * diff[nonzero] / np.abs(ref_values[nonzero])))
            if np.any(nonzero) else 0.0
        )
    return np.array(nus), np.array(absolute), np.array(relative)


def plot_case(case_dir: Path, output: Path | None = None) -> None:
    case_dir = case_dir.expanduser().resolve()
    validation_root = case_dir.parents[2]
    reference_path = validation_root / "references" / case_dir.name / "rad.nc"
    candidate_path = case_dir / "rad.nc"
    for path in (reference_path, candidate_path):
        if not path.is_file():
            raise FileNotFoundError(path)
    accuracy_path = case_dir / "accuracy.json"
    accuracy = json.loads(accuracy_path.read_text()) if accuracy_path.is_file() else {}
    relative_tolerance = float(accuracy.get("threshold_percent", 0.01))
    science_abs_tolerance = float(accuracy.get("bt_abs_tolerance", 1e-10))
    tau_abs_tolerance = float(accuracy.get("tau_abs_tolerance", 1e-12))

    with nc.Dataset(reference_path) as reference, nc.Dataset(candidate_path) as candidate:
        science_prefixes = {
            name.split("_", 1)[0] + "_" for name in reference.variables
            if name.startswith(("bt_", "rad_"))
        }
        if len(science_prefixes) != 1:
            raise ValueError(f"expected one science-variable family, found {sorted(science_prefixes)}")
        science_prefix = science_prefixes.pop()
        science = error_series(reference, candidate, science_prefix)
        tau = error_series(reference, candidate, "tau_")

    fig, axes = plt.subplots(2, 2, figsize=(13, 8), sharex="col", constrained_layout=True)
    configurations = (
        (axes[0, 0], science[0], science[1], science_abs_tolerance, "Science absolute difference"),
        (axes[1, 0], science[0], science[2], relative_tolerance, "Science relative difference (%)"),
        (axes[0, 1], tau[0], tau[1], tau_abs_tolerance, "Transmittance absolute difference"),
        (axes[1, 1], tau[0], tau[2], relative_tolerance, "Transmittance relative difference (%)"),
    )
    for ax, nus, values, tolerance, title in configurations:
        ax.plot(nus, values, "o-", color="steelblue", markersize=3, linewidth=1)
        ax.axhline(tolerance, color="darkorange", linestyle="--", linewidth=1.2, label="Tolerance")
        ax.set_title(title)
        ax.set_ylabel("Difference")
        ax.grid(alpha=0.3)
        ax.legend()
        ymax = max(float(np.nanmax(values)) if values.size else 0.0, tolerance)
        ax.set_ylim(0, max(1e-18, 1.15 * ymax))
    axes[1, 0].set_xlabel("Wavenumber (cm$^{-1}$)")
    axes[1, 1].set_xlabel("Wavenumber (cm$^{-1}$)")
    quantity = "Radiance" if science_prefix == "rad_" else "Brightness temperature"
    fig.suptitle(f"{case_dir.name}: {quantity} and transmittance errors")
    selected_output = output or case_dir / "channel_errors.png"
    selected_output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(selected_output, dpi=150)
    plt.close(fig)
    print(f"Saved: {selected_output}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case_dir", type=Path)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    plot_case(args.case_dir, args.output)
