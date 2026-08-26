#!/usr/bin/env python3
"""Plot channel-wise validation errors in geometry-specific spectral panels."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RUNS = ROOT / "runs"
REFERENCES = ROOT / "references"
GEOMETRIES = ("limb", "nadir", "zenith")


def latest_run() -> Path:
    candidates = sorted(path for path in RUNS.glob("validation_*") if (path / "summary.tsv").is_file())
    if not candidates:
        raise FileNotFoundError(f"no validation summary found below {RUNS}")
    return candidates[-1]


def resolve_run_dir(value: str | None) -> Path:
    configured = value or os.environ.get("VALIDATION_RUN_DIR")
    return Path(configured).expanduser().resolve() if configured else latest_run().resolve()


def case_tolerance(run_dir: Path, case_name: str, fallback: float) -> float:
    path = run_dir / case_name / "accuracy.json"
    if path.is_file():
        try:
            return float(json.loads(path.read_text()).get("threshold_percent", fallback))
        except (OSError, ValueError, json.JSONDecodeError):
            pass
    return fallback


def error_statistics(reference, candidate, prefix: str, relative: bool):
    ref_names = sorted(
        (name for name in reference.variables if name.startswith(prefix)),
        key=lambda name: float(name.split("_", 1)[1]),
    )
    test_names = sorted(
        (name for name in candidate.variables if name.startswith(prefix)),
        key=lambda name: float(name.split("_", 1)[1]),
    )
    if ref_names != test_names:
        raise ValueError(f"reference and candidate {prefix} variable sets differ")
    nus = []
    statistics = []
    for name in ref_names:
        ref = np.ma.filled(reference.variables[name][:], np.nan).astype(float)
        test = np.ma.filled(candidate.variables[name][:], np.nan).astype(float)
        if ref.shape != test.shape:
            raise ValueError(f"shape mismatch for {name}: {ref.shape} != {test.shape}")
        ref_finite = np.isfinite(ref)
        test_finite = np.isfinite(test)
        if not np.array_equal(ref_finite, test_finite):
            raise ValueError(f"finite/NaN mask mismatch for {name}")
        usable = ref_finite & test_finite
        if relative:
            usable &= ref != 0
            errors = 100.0 * (test[usable] - ref[usable]) / ref[usable]
        else:
            errors = test[usable] - ref[usable]
        nus.append(float(name.split("_", 1)[1]))
        if errors.size:
            statistics.append((
                float(np.mean(errors)),
                float(np.std(errors)),
                float(np.min(errors)),
                float(np.max(errors)),
            ))
        else:
            statistics.append((0.0, 0.0, 0.0, 0.0))
    return np.array(nus), np.array(statistics)


def merge_statistics(
    target: dict[float, np.ndarray], nus: np.ndarray, statistics: np.ndarray
) -> None:
    for nu, values in zip(nus, statistics):
        key = float(nu)
        if key in target and not np.allclose(target[key], values, equal_nan=True):
            raise ValueError(f"overlapping validation cases disagree at {key:.4f} cm^-1")
        target.setdefault(key, values)


def collect_geometry(run_dir: Path, subset: pd.DataFrame, fallback_tolerance: float):
    science_re: dict[float, np.ndarray] = {}
    tau_re: dict[float, np.ndarray] = {}
    science_ae: dict[float, np.ndarray] = {}
    tau_ae: dict[float, np.ndarray] = {}
    tolerances: dict[float, float] = {}
    for case_name in subset["case_name"]:
        reference_path = REFERENCES / case_name / "rad.nc"
        candidate_path = run_dir / case_name / "rad.nc"
        for required in (reference_path, candidate_path):
            if not required.is_file():
                raise FileNotFoundError(required)
        tolerance = case_tolerance(run_dir, case_name, fallback_tolerance)
        with nc.Dataset(reference_path) as reference, nc.Dataset(candidate_path) as candidate:
            science_prefixes = {
                name.split("_", 1)[0] + "_" for name in reference.variables
                if name.startswith(("bt_", "rad_"))
            }
            if len(science_prefixes) != 1:
                raise ValueError(
                    f"expected one science-variable family in {reference_path}, "
                    f"found {sorted(science_prefixes)}"
                )
            science_prefix = science_prefixes.pop()
            science_nus, science_re_statistics = error_statistics(
                reference, candidate, science_prefix, relative=True
            )
            _, science_ae_statistics = error_statistics(
                reference, candidate, science_prefix, relative=False
            )
            tau_nus, tau_re_statistics = error_statistics(
                reference, candidate, "tau_", relative=True
            )
            _, tau_ae_statistics = error_statistics(
                reference, candidate, "tau_", relative=False
            )
        merge_statistics(science_re, science_nus, science_re_statistics)
        merge_statistics(tau_re, tau_nus, tau_re_statistics)
        merge_statistics(science_ae, science_nus, science_ae_statistics)
        merge_statistics(tau_ae, tau_nus, tau_ae_statistics)
        for nu in np.concatenate((science_nus, tau_nus)):
            value = float(nu)
            tolerances[value] = min(tolerances.get(value, tolerance), tolerance)
    nus = np.array(sorted(set(science_re) | set(tau_re)))
    return (
        nus,
        np.stack([science_re[float(nu)] for nu in nus]),
        np.stack([tau_re[float(nu)] for nu in nus]),
        np.stack([science_ae[float(nu)] for nu in nus]),
        np.stack([tau_ae[float(nu)] for nu in nus]),
        np.array([tolerances[float(nu)] for nu in nus]),
    )


def plot_statistics(ax, nus, statistics, tolerance, title: str, suffix: str) -> None:
    styles = (
        (f"Mean{suffix}", "steelblue", "-"),
        (f"SD{suffix}", "mediumpurple", "-"),
        (f"Min{suffix}", "seagreen", "-"),
        (f"Max{suffix}", "crimson", "-"),
    )
    for index, (label, color, linestyle) in enumerate(styles):
        ax.plot(nus, statistics[:, index], color=color, linestyle=linestyle,
                marker="o", markersize=2, linewidth=0.9, label=label)
    if tolerance is not None:
        ax.plot(nus, tolerance, "--", color="darkorange", linewidth=1.1,
                label="± tolerance")
        ax.plot(nus, -tolerance, "--", color="darkorange", linewidth=1.1)
    ax.axhline(0, color="0.45", linewidth=0.8, zorder=0)
    ax.set_title(title)
    ax.grid(alpha=0.3)
    ax.legend(loc="upper right", ncol=2, fontsize=8)
    values = statistics.ravel() if tolerance is None else np.concatenate((statistics.ravel(), tolerance))
    extent = float(np.nanmax(np.abs(values)))
    span = max(1e-12, extent)
    ax.set_ylim(-1.2 * span, 1.2 * span)


def plot_summary(
    run_dir: Path,
    combined_output: Path,
    relative_output: Path,
    absolute_output: Path,
    fallback_tolerance: float,
) -> None:
    summary = run_dir / "summary.tsv"
    if not summary.is_file():
        raise FileNotFoundError(f"summary not found: {summary}")
    frame = pd.read_csv(summary, sep="\t")
    required = {"case_name", "geometry", "status"}
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"summary lacks columns: {', '.join(sorted(missing))}")

    figures = {
        "combined": plt.subplots(3, 2, figsize=(15, 11), sharex=True, constrained_layout=True),
        "RE": plt.subplots(3, 2, figsize=(15, 11), sharex=True, constrained_layout=True),
        "AE": plt.subplots(3, 2, figsize=(15, 11), sharex=True, constrained_layout=True),
    }
    for row, geometry in enumerate(GEOMETRIES):
        subset = frame[frame["geometry"] == geometry].sort_values("case_name")
        if subset.empty:
            for _, axes in figures.values():
                for ax in axes[row]:
                    ax.text(0.5, 0.5, "No cases", ha="center", va="center", transform=ax.transAxes)
            continue
        nus, science_re, tau_re, science_ae, tau_ae, tolerance = collect_geometry(
            run_dir, subset, fallback_tolerance
        )
        science_label = "Radiance" if geometry == "limb" else "Brightness temperature"
        combined_axes = figures["combined"][1]
        re_axes = figures["RE"][1]
        ae_axes = figures["AE"][1]
        plot_statistics(
            re_axes[row, 0], nus, science_re, tolerance,
            f"{geometry.capitalize()}: {science_label} ({len(nus)} channels)", "RE",
        )
        plot_statistics(
            re_axes[row, 1], nus, tau_re, tolerance,
            f"{geometry.capitalize()}: Transmittance τ ({len(nus)} channels)", "RE",
        )
        plot_statistics(
            ae_axes[row, 0], nus, science_ae, None,
            f"{geometry.capitalize()}: {science_label} ({len(nus)} channels)", "AE",
        )
        plot_statistics(
            ae_axes[row, 1], nus, tau_ae, None,
            f"{geometry.capitalize()}: Transmittance τ ({len(nus)} channels)", "AE",
        )
        if geometry == "limb":
            plot_statistics(
                combined_axes[row, 0], nus, science_re, tolerance,
                f"{geometry.capitalize()}: {science_label} RE ({len(nus)} channels)", "RE",
            )
            combined_axes[row, 0].set_ylabel("Relative error (%)")
        else:
            plot_statistics(
                combined_axes[row, 0], nus, science_ae, None,
                f"{geometry.capitalize()}: {science_label} AE ({len(nus)} channels)", "AE",
            )
            combined_axes[row, 0].set_ylabel("Absolute error (K)")
        plot_statistics(
            combined_axes[row, 1], nus, tau_ae, None,
            f"{geometry.capitalize()}: Transmittance τ AE ({len(nus)} channels)", "AE",
        )
        combined_axes[row, 1].set_ylabel("Absolute error")
        re_axes[row, 0].set_ylabel("Relative error (%)")
        science_unit = "W m$^{-2}$ sr$^{-1}$ (cm$^{-1}$)$^{-1}$" if geometry == "limb" else "K"
        ae_axes[row, 0].set_ylabel(f"Absolute error ({science_unit})")
        ae_axes[row, 1].set_ylabel("Absolute error")
    profile = ", ".join(sorted(frame["profile"].astype(str).unique())) if "profile" in frame else "validation"
    for suffix, (fig, axes), output in (
        ("selected metrics", figures["combined"], combined_output),
        ("relative error", figures["RE"], relative_output),
        ("absolute error", figures["AE"], absolute_output),
    ):
        axes[-1, 0].set_xlabel("Wavenumber (cm$^{-1}$)")
        axes[-1, 1].set_xlabel("Wavenumber (cm$^{-1}$)")
        fig.suptitle(f"Validation summary ({suffix}): {run_dir.name} ({profile})", fontsize=14)
        output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output, dpi=150)
        plt.close(fig)
        print(f"Saved: {output}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_dir", nargs="?", help="Defaults to VALIDATION_RUN_DIR or the latest run.")
    parser.add_argument("--output", type=Path, help="Relative-error plot output.")
    parser.add_argument("--absolute-output", type=Path, help="Absolute-error plot output.")
    parser.add_argument("--combined-output", type=Path, help="Combined summary plot output.")
    parser.add_argument("--default-tolerance", type=float,
                        default=float(os.environ.get("VALIDATION_DEFAULT_TOLERANCE_PERCENT", "0.01")))
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    selected_run = resolve_run_dir(args.run_dir)
    plot_summary(
        selected_run,
        args.combined_output or selected_run / "summary_plot.png",
        args.output or selected_run / "summary_plot_re.png",
        args.absolute_output or selected_run / "summary_plot_ae.png",
        args.default_tolerance,
    )
