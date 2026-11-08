#!/usr/bin/env python3
"""Summarize the synthetic retrieval and plot profiles and radiance residuals."""

import csv
import json
import math
from pathlib import Path
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, TwoSlopeNorm
import numpy as np

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
CASE = WORK / "case"
# Keep the control-file order for ASCII radiance columns; matrix plots sort it.
CHANNELS = (791, 792, 793, 810, 811, 812, 778, 779, 780,
            787, 788, 863, 864, 867, 869)
# Key, JURASSIC state name, zero-based atmosphere column, plot scale/label,
# and derivative unit for each retrieved quantity.
QUANTITIES = (("temperature", "TEMPERATURE", 5, 1, "Temperature (K)", "K"),
              ("o3", "O3", 8, 1e6, "O3 (ppmv)", "ppmv"),
              ("hno3", "HNO3", 9, 1e9, "HNO3 (ppbv)", "ppbv"))
# Representative configured channels for signed kernel slices (cm^-1).
SLICE_CHANNELS = {"temperature": 792, "o3": 788, "hno3": 867}


def rows(path):
    """Read numeric ASCII rows, skipping the JURASSIC column descriptions.

    JURASSIC headers number columns from 1; Python lists index from 0.
    Atmosphere columns 2/6/9/10 (Python 1/5/8/9) hold z/T/O3/HNO3.
    Observation column 5 (Python 4) is view altitude, and radiances
    start at column 11 (Python 10).
    """
    return [[float(value) for value in line.split()]
            for line in path.read_text().splitlines()
            if line and not line.startswith("#")]


def rms(values):
    """Return root mean square of equally weighted values."""
    return math.sqrt(sum(value * value for value in values) / len(values))


def timer(log, name):
    """Extract mean elapsed seconds from a JURASSIC timer line."""
    match = re.search(rf"^TIMER_{name} = .*?mean=\s*([\d.eE+-]+) s",
                      log, re.MULTILINE)
    return float(match.group(1)) if match else None


def read_matrix(path, row_space, col_space, shape):
    """Reconstruct a dense MATRIXFMT=1 matrix from its ASCII element records.

    Each text row is one matrix element: columns 1--6 describe its matrix
    row, 7--12 its matrix column, and 13 its value. Explicit integer row
    and column indices reconstruct the array regardless of textual order.
    Retained state metadata identify quantity and altitude; measurement
    metadata identify spectral channel and tangent/view altitude for checks.
    """
    lines = path.read_text().splitlines()
    header = "\n".join(line for line in lines if line.startswith("#"))
    for number, space, side in ((1, row_space, "Row"), (7, col_space, "Col")):
        label = f"# ${number} = {side}: index ({space} space)"
        if label not in header:
            raise ValueError(f"Unexpected matrix header in {path}: {label}")
    if "# $13 = Matrix element" not in header:
        raise ValueError(f"Missing value header in {path}")
    matrix = np.full(shape, np.nan)
    row_meta, col_meta = {}, {}
    for line in lines:
        if not line or line.startswith("#"):
            continue
        fields = line.split()
        if len(fields) != 13:
            raise ValueError(f"Expected 13 fields in {path}: {line}")
        # Columns 1 and 7 are zero-based matrix indices, unlike the header's
        # one-based numbering of ASCII fields; reject gaps or duplicates.
        i, j = int(fields[0]), int(fields[6])
        if not (0 <= i < shape[0] and 0 <= j < shape[1]) or not np.isnan(matrix[i, j]):
            raise ValueError(f"Invalid or repeated matrix index {(i, j)} in {path}")
        # Repeated row/column coordinates must carry identical metadata.
        for metadata, index, values in ((row_meta, i, tuple(fields[1:6])),
                                        (col_meta, j, tuple(fields[7:12]))):
            if index in metadata and metadata[index] != values:
                raise ValueError(f"Inconsistent matrix metadata at index {index} in {path}")
            metadata[index] = values
        matrix[i, j] = float(fields[12])
    if not np.isfinite(matrix).all() or len(row_meta) != shape[0] or len(col_meta) != shape[1]:
        raise ValueError(f"Incomplete or non-finite matrix: {path}")
    return matrix, [row_meta[i] for i in range(shape[0])], [col_meta[j] for j in range(shape[1])]


def matrix_diagnostics(prior, final, selected, metrics):
    """Diagnose spectral sensitivity, vertical response, and formal uncertainty."""
    # Read the finite-difference kernel, averaging kernel, and a priori and
    # posterior covariance matrices written by this retrieval.
    n_state = 3 * len(selected)
    kernel, measurements, states = read_matrix(
        CASE / "matrix_kernel.tab", "measurement", "state", (16 * len(CHANNELS), n_state))
    avk, avk_rows, avk_cols = read_matrix(
        CASE / "matrix_avk.tab", "state", "state", (n_state, n_state))
    cov_apr, apr_rows, apr_cols = read_matrix(
        CASE / "matrix_cov_apr.tab", "state", "state", (n_state, n_state))
    cov_ret, ret_rows, ret_cols = read_matrix(
        CASE / "matrix_cov_ret.tab", "state", "state", (n_state, n_state))
    # A common state-vector ordering is required to compare matrix blocks.
    if not (states == avk_rows == avk_cols == apr_rows == apr_cols == ret_rows == ret_cols):
        raise ValueError("State metadata differ between retrieval matrices")
    # Check every measurement row against the configured spectral channels
    # and 16 limb altitudes before grouping rows for sensitivity diagnostics.
    if sorted(float(m[0]) for m in measurements) != sorted(
            float(nu) for nu in CHANNELS for _ in range(16)):
        raise ValueError("Kernel channels differ from the control file")
    if sorted(float(m[2]) for m in measurements) != sorted(
            float(z) for z in range(15, 46, 2) for _ in CHANNELS):
        raise ValueError("Kernel view altitudes differ from limb geometry")

    fig_jac, jac_axes = plt.subplots(1, 3, figsize=(13, 5), sharey=True)
    fig_slice, slice_axes = plt.subplots(1, 3, figsize=(12, 5), sharey=True)
    fig_avk, avk_axes = plt.subplots(1, 3, figsize=(12, 5), sharey=True)
    fig_unc, unc_axes = plt.subplots(1, 3, figsize=(12, 5), sharey=True)
    channels = sorted(CHANNELS)
    for ax_jac, ax_slice, ax_avk, ax_unc, (key, name, column, scale, _, unit) in zip(
            jac_axes, slice_axes, avk_axes, unc_axes, QUANTITIES):
        # State columns carry quantity labels and altitude coordinates. Select
        # one 16-level block and verify it matches the retrieved profile grid.
        state_ids = [j for j, meta in enumerate(states) if meta[0] == name]
        altitudes = np.array([float(states[j][2]) for j in state_ids])
        if len(state_ids) != len(selected) or not np.array_equal(
                altitudes, [prior[i][1] for i in selected]):
            raise ValueError(f"Unexpected state grid for {name}")
        edges = np.r_[altitudes - 1, altitudes[-1] + 1]

        # Kernel rows are channel/limb measurements; columns are state levels.
        # A channel occurs once at each of 16 tangent altitudes. For a fixed
        # channel and state level, RMS over those rows is a positive sensitivity
        # magnitude, so it loses derivative sign. Dividing by the profile plot
        # scale converts per-unit-VMR derivatives to per-ppmv or per-ppbv.
        sensitivity = np.array([
            np.sqrt(np.mean(kernel[np.ix_(
                [i for i, meta in enumerate(measurements) if float(meta[0]) == nu],
                state_ids)] ** 2, axis=0)) / scale for nu in channels]).T
        positive = sensitivity[sensitivity > 0]
        if len(positive) == 0:
            raise ValueError(f"Zero finite-difference sensitivity for {name}")
        maximum = float(np.max(positive))
        image = ax_jac.pcolormesh(np.arange(len(channels) + 1), edges,
                                  np.ma.masked_less_equal(sensitivity, 0), cmap="viridis",
                                  norm=LogNorm(vmin=max(float(np.min(positive)), maximum * 1e-4),
                                               vmax=maximum))
        ax_jac.set(xlabel="Wavenumber (cm⁻¹)", ylim=(15, 45),
                   xticks=np.arange(len(channels)) + 0.5, title=name)
        ax_jac.set_xticklabels(channels, rotation=90)
        fig_jac.colorbar(image, ax=ax_jac, label=f"W/(m² sr cm⁻¹) per {unit}")

        # Show one signed channel slice directly from measurement rows and
        # state columns. Transpose to put state altitude on the vertical axis.
        nu = SLICE_CHANNELS[key]
        ray_ids = sorted((i for i, meta in enumerate(measurements)
                          if float(meta[0]) == nu),
                         key=lambda i: float(measurements[i][2]))
        view_altitudes = np.array([float(measurements[i][2]) for i in ray_ids])
        if not np.array_equal(view_altitudes, altitudes):
            raise ValueError(f"Unexpected limb grid at {nu} cm^-1")
        view_edges = np.r_[view_altitudes - 1, view_altitudes[-1] + 1]
        signed = (kernel[np.ix_(ray_ids, state_ids)] / scale).T
        signed_bound = float(np.max(np.abs(signed)))
        if signed_bound <= 0:
            raise ValueError(f"Zero signed kernel slice for {name} at {nu} cm^-1")
        image = ax_slice.pcolormesh(view_edges, edges, signed, cmap="RdBu_r",
                                    norm=TwoSlopeNorm(vmin=-signed_bound, vcenter=0,
                                                      vmax=signed_bound))
        ax_slice.set(xlabel="Tangent altitude (km)", title=f"{name} · {nu} cm⁻¹",
                     xlim=(15, 45), ylim=(15, 45))
        fig_slice.colorbar(image, ax=ax_slice, label=f"W/(m² sr cm⁻¹) per {unit}")

        # The within-quantity block shows vertical response in comparable
        # units; cross-quantity blocks remain in the full matrix but are not
        # shown here. Its trace is DOFS, and each row sum is measurement
        # response, matching JURASSIC's analyze_avk_quantity convention.
        block = avk[np.ix_(state_ids, state_ids)]
        metrics[f"{key}_dofs"] = float(np.trace(block))
        metrics[f"{key}_measurement_response_mean"] = float(np.mean(np.sum(block, axis=1)))
        bound = float(np.max(np.abs(block)))
        image = ax_avk.pcolormesh(edges, edges, block, cmap="RdBu_r",
                                  norm=TwoSlopeNorm(vmin=-bound, vcenter=0, vmax=bound))
        ax_avk.set(xlabel="State altitude (km)", xlim=(15, 45),
                   ylim=(15, 45), title=name)
        fig_avk.colorbar(image, ax=ax_avk, label="Averaging-kernel value")

        # Covariance diagonals are variances; square roots give formal 1-sigma
        # uncertainty. Temperature stays in K. Gas percentages use the
        # corresponding a priori or retrieved VMR as their denominator.
        apr_var = np.diag(cov_apr)[state_ids]
        ret_var = np.diag(cov_ret)[state_ids]
        if np.any(apr_var <= 0) or np.any(ret_var <= 0):
            raise ValueError(f"Nonpositive covariance diagonal for {name}")
        prior_sigma, final_sigma = np.sqrt(apr_var), np.sqrt(ret_var)
        if key != "temperature":
            prior_vmr = np.array([prior[i][column] for i in selected])
            final_vmr = np.array([final[i][column] for i in selected])
            if np.any(prior_vmr <= 0) or np.any(final_vmr <= 0):
                raise ValueError(f"Nonpositive gas VMR for {name}")
            prior_sigma = 100 * prior_sigma / prior_vmr
            final_sigma = 100 * final_sigma / final_vmr
        ax_unc.plot(prior_sigma, altitudes, "o-", markersize=3.5, label="A priori")
        ax_unc.plot(final_sigma, altitudes, "o-", markersize=3.5, label="Posterior")
        ax_unc.set(xlabel="1σ uncertainty (K)" if key == "temperature"
                   else "Relative 1σ uncertainty (%)", title=name, ylim=(15, 45))
        ax_unc.grid(alpha=0.3)
        # Median fractional reduction compares posterior and prior standard
        # deviations in the original state units, before gas normalization.
        metrics[f"{key}_median_uncertainty_reduction_percent"] = float(
            100 * np.median(1 - np.sqrt(ret_var / apr_var)))

    unc_axes[0].legend()
    for fig, axes, filename, ylabel, title in (
            (fig_jac, jac_axes, "jacobians.png", "Retrieval altitude (km)",
             "RMS finite-difference sensitivity over limb observations"),
            (fig_slice, slice_axes, "kernel_slices.png", "Retrieval altitude (km)",
             "Signed finite-difference kernel at representative channels"),
            (fig_avk, avk_axes, "averaging_kernels.png", "Retrieval altitude (km)",
             "Averaging kernels (within-quantity blocks)"),
            (fig_unc, unc_axes, "uncertainties.png", "Altitude (km)", None)):
        axes[0].set_ylabel(ylabel)
        if title:
            fig.suptitle(title)
        fig.tight_layout()
        fig.savefig(WORK / filename, dpi=160)
        plt.close(fig)


def main():
    """Calculate profile and radiance errors, then write plots and matrices."""
    # Read the unperturbed a priori, synthetic truth, retrieved atmosphere,
    # and model/observation radiances from the completed example.
    prior = rows(CASE / "atm_apr.tab")
    truth = rows(CASE / "atm_truth.tab")
    final = rows(CASE / "atm_final.tab")
    initial_radiance = rows(CASE / "obs_initial.tab")
    measured = rows(CASE / "obs_meas.tab")
    final_radiance = rows(CASE / "obs_final.tab")
    assert len(prior) == len(truth) == len(final)
    assert len(initial_radiance) == len(measured) == len(final_radiance) == 16
    # Restrict profile diagnostics to the 15--45 km retrieval state grid.
    selected = [i for i, row in enumerate(truth) if 15 <= row[1] <= 45]
    metrics = {}
    # Profile RMSE covers those 16 levels. Gas errors are levelwise relative
    # percent errors against synthetic truth, then combined by RMS.
    for label, _, column, _, _, _ in QUANTITIES:
        if label == "temperature":
            error = lambda value, target: value - target
            suffix = "rmse_K"
        else:
            error = lambda value, target: 100 * (value / target - 1)
            suffix = "relative_rmse_percent"
        metrics[f"{label}_prior_{suffix}"] = rms(
            [error(prior[i][column], truth[i][column]) for i in selected])
        metrics[f"{label}_retrieved_{suffix}"] = rms(
            [error(final[i][column], truth[i][column]) for i in selected])

    # Initial and final residuals both use observation - model. Observation
    # radiances start at JURASSIC ASCII column 11 in control-file order.
    residuals = []
    for initial, observation, result in zip(initial_radiance, measured, final_radiance):
        assert initial[4] == observation[4] == result[4]
        for channel, nu in enumerate(CHANNELS):
            residuals.append((observation[4], nu, observation[10 + channel] - initial[10 + channel],
                              observation[10 + channel] - result[10 + channel]))
    metrics["initial_radiance_residual_rms"] = rms([r[2] for r in residuals])
    metrics["final_radiance_residual_rms"] = rms([r[3] for r in residuals])
    # JURASSIC reports chi^2/m (cost per measurement) for each iteration.
    # Preserve its iteration history for the convergence plot.
    log = (WORK / "retrieval.log").read_text()
    iterations = re.findall(r"it=\s*(\d+) / chi\^2/m=\s*([\d.eE+-]+)", log)
    metrics["iterations"] = int(iterations[-1][0]) if iterations else None
    metrics["converged"] = metrics["iterations"] is not None and metrics["iterations"] < 30
    metrics["final_cost_per_measurement"] = float(iterations[-1][1]) if iterations else None
    # Wall time ends after retrieval, before plotting. RET_ITERATE includes
    # both optimization and later Jacobian recomputations, not kernel time alone.
    metrics["total_runtime_s"] = float((WORK / "total_runtime_s.txt").read_text())
    metrics["retrieval_kernel_initial_s"] = timer(log, "RET_KERNEL_INIT")
    metrics["retrieval_iterate_s"] = timer(log, "RET_ITERATE")
    metrics["retrieval_setup_s"] = timer(log, "RET_SETUP")
    metrics["retrieval_diagnostics_s"] = timer(log, "RET_DIAGNOSTICS")
    # Keep retrieved open markers above the truth line so their close
    # agreement remains visible without moving either profile.
    fig, axes = plt.subplots(1, 3, figsize=(12, 5), sharey=True)
    for ax, (_, _, column, scale, title, _) in zip(axes, QUANTITIES):
        for name, table, style in (
                ("Truth", truth, {"linestyle": "-", "linewidth": 2, "zorder": 1}),
                ("A priori", prior, {"linestyle": "--", "marker": "s",
                                    "markersize": 3.5, "zorder": 2}),
                ("Retrieved", final, {"linestyle": ":", "marker": "o",
                                     "markersize": 6, "markerfacecolor": "white",
                                     "markeredgewidth": 1.4, "zorder": 3})):
            ax.plot([table[i][column] * scale for i in selected],
                    [table[i][1] for i in selected], label=name, **style)
        ax.set_xlabel(title)
        ax.grid(alpha=0.3)
        ax.set_ylim(15, 45)
    axes[0].set_ylabel("Altitude (km)")
    axes[0].legend()
    fig.tight_layout()
    fig.savefig(WORK / "profiles.png", dpi=160)
    plt.close(fig)

    # Plot residual RMS across channels by tangent altitude and the cost
    # history by iteration, with logarithmic residual/cost axes.
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    heights = sorted({r[0] for r in residuals})
    for label, index in (("Initial", 2), ("Final", 3)):
        axes[0].semilogx([rms([r[index] for r in residuals if r[0] == z])
                         for z in heights], heights, "o-", label=label)
    axes[0].set(xlabel="Radiance residual RMS across all channels\n[W/(m² sr cm⁻¹)]",
                ylabel="Tangent altitude (km)", ylim=(15, 45))
    axes[0].legend()
    if iterations:
        axes[1].semilogy([int(it) for it, _ in iterations], [float(cost) for _, cost in iterations], "o-")
    axes[1].set(xlabel="Iteration", ylabel="Cost / measurement")
    for ax in axes:
        ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(WORK / "residuals_convergence.png", dpi=160)
    plt.close(fig)

    # Generate finite-difference, averaging-kernel, and covariance diagnostics.
    matrix_diagnostics(prior, final, selected, metrics)
    # Store machine-readable metrics and every signed radiance residual.
    (WORK / "metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    with (WORK / "residuals.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("view_altitude_km", "wavenumber_cm-1", "initial_residual", "final_residual"))
        writer.writerows(residuals)
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
