#!/usr/bin/env python3
"""Create geometry-appropriate accuracy plots from a completed validation run."""

import argparse
import csv
import math
from pathlib import Path

HEIGHTS = (5, 10, 20, 50)
METHODS = (("ega", "EGA", "#0072B2"), ("cga", "CGA", "#E69F00"))
C1 = 1.19104259e-8
C2 = 1.43877506


def load_spectrum(path, geometry):
    """Load one geometry keyed by ray number and channel wavenumber."""
    values = {}
    if not path.is_file():
        raise FileNotFoundError(f"missing spectrum: {path}")
    with path.open(newline="") as stream:
        for row in csv.DictReader(stream):
            if row["geometry"] == geometry:
                values[(int(row["ray"]), float(row["nu_cm-1"]))] = float(
                    row["radiance_W_m-2_sr-1_cm"])
    if not values:
        raise ValueError(f"no {geometry} spectra in {path}")
    return values


def read_rows(root, geometry, method):
    """Align one JURASSIC method with RFM and calculate channel errors."""
    reference = load_spectrum(root / "rfm_reference" / "spectra.csv", geometry)
    candidate = load_spectrum(root / f"test_{method}" / "spectra.csv", geometry)
    if candidate.keys() != reference.keys():
        raise ValueError(f"{method}/RFM grid mismatch for {geometry}")
    rows = []
    for (ray, nu), rfm in sorted(reference.items()):
        value = candidate[ray, nu]
        difference = value - rfm
        relative = math.nan if rfm == 0 else 100 * difference / rfm
        rows.append({"nu_cm-1": str(nu), "ray": str(ray),
                     "jurassic_radiance": str(value), "rfm_radiance": str(rfm),
                     "signed_difference": str(difference),
                     "absolute_difference": str(abs(difference)),
                     "relative_difference_percent": str(relative)})
    return rows


def brightness_temperature(radiance, wavenumber):
    """JURASSIC BRIGHT(rad, nu), applied after channel averaging."""
    if radiance <= 0:
        return math.nan
    return C2 * wavenumber / math.log1p(C1 * wavenumber ** 3 / radiance)


def percentile(values, fraction):
    """Return a linearly interpolated percentile from a nonempty sequence."""
    values = sorted(values)
    position = fraction * (len(values) - 1)
    lower = int(position)
    weight = position - lower
    return values[lower] * (1 - weight) + values[min(lower + 1, len(values) - 1)] * weight


def statistics(values):
    """Summarize finite absolute errors without discarding large values."""
    absolute = [abs(value) for value in values if math.isfinite(value)]
    return {
        "count": len(absolute),
        "rms": math.sqrt(sum(value * value for value in absolute) / len(absolute)),
        "median": percentile(absolute, 0.5),
        "p95": percentile(absolute, 0.95),
        "maximum": max(absolute),
    }


def limb_data(root, method):
    rows = read_rows(root, "limb", method)
    return [[row for row in rows if int(row["ray"]) == ray] for ray in range(len(HEIGHTS))]


def limb_statistics(rows):
    """Summarize relative errors and retain context for the largest value."""
    result = statistics(float(row["relative_difference_percent"]) for row in rows)
    maximum = max(rows, key=lambda row: abs(float(row["relative_difference_percent"])))
    result.update({
        "maximum_channel_cm-1": float(maximum["nu_cm-1"]),
        "reference_at_maximum": float(maximum["rfm_radiance"]),
        "absolute_difference_at_maximum": float(maximum["absolute_difference"]),
    })
    return result


def bt_data(root, geometry, method):
    rows = read_rows(root, geometry, method)
    result = []
    for row in rows:
        nu = float(row["nu_cm-1"])
        jurassic = brightness_temperature(float(row["jurassic_radiance"]), nu)
        rfm = brightness_temperature(float(row["rfm_radiance"]), nu)
        result.append((nu, jurassic, rfm, jurassic - rfm))
    return result


def plot_limb_radiance(root, output, plt):
    data = {method: limb_data(root, method) for method, _, _ in METHODS}
    fig, axes = plt.subplots(4, 1, figsize=(11, 9), sharex=True, constrained_layout=True)
    for ray, (ax, height) in enumerate(zip(axes, HEIGHTS)):
        reference = data["ega"][ray]
        nu = [float(row["nu_cm-1"]) for row in reference]
        ax.plot(nu, [float(row["rfm_radiance"]) for row in reference], color="0.25",
                lw=0.75, label="RFM")
        for method, label, color in METHODS:
            rows = data[method][ray]
            ax.plot(nu, [float(row["jurassic_radiance"]) for row in rows], color=color,
                    lw=0.6, label=f"JURASSIC {label}")
        ax.set_yscale("log")
        ax.set_ylabel("Radiance\n[W m$^{-2}$ sr$^{-1}$ cm]")
        ax.set_title(f"Limb: {height} km geometric tangent height", loc="left")
        ax.grid(alpha=0.18, lw=0.5)
    axes[0].legend(frameon=False, ncol=4, fontsize=8)
    axes[-1].set_xlabel("Wavenumber [cm$^{-1}$]")
    fig.suptitle("Limb radiance spectra")
    fig.savefig(output, dpi=220)
    plt.close(fig)


def plot_limb_relative(root, output, plt):
    data = {method: limb_data(root, method) for method, _, _ in METHODS}
    fig, axes = plt.subplots(4, 1, figsize=(11, 9), sharex=True, constrained_layout=True)
    for ray, (ax, height) in enumerate(zip(axes, HEIGHTS)):
        for method, label, color in METHODS:
            rows = data[method][ray]
            nu = [float(row["nu_cm-1"]) for row in rows]
            relative = [float(row["relative_difference_percent"]) for row in rows]
            ax.plot(nu, relative, color=color, lw=0.6, label=f"{label} − RFM")
        ax.axhline(0, color="0.25", lw=0.5)
        ax.set_ylabel("Difference [%]")
        ax.set_title(f"Limb: {height} km geometric tangent height", loc="left")
        ax.grid(alpha=0.18, lw=0.5)
    axes[0].legend(frameon=False, ncol=2)
    axes[-1].set_xlabel("Wavenumber [cm$^{-1}$]")
    fig.suptitle("Relative limb radiance differences")
    fig.savefig(output, dpi=220)
    plt.close(fig)


def plot_bt(root, geometry, output, plt):
    method_data = {method: bt_data(root, geometry, method) for method, _, _ in METHODS}
    reference = method_data["ega"]
    nu = [row[0] for row in reference]
    fig, axes = plt.subplots(2, 1, figsize=(11, 6.5), sharex=True, constrained_layout=True)
    axes[0].plot(nu, [row[2] for row in reference], color="0.25", lw=0.75, label="RFM")
    for method, label, color in METHODS:
        data = method_data[method]
        axes[0].plot(nu, [row[1] for row in data], color=color, lw=0.6,
                     label=f"JURASSIC {label}")
        axes[1].plot(nu, [row[3] for row in data], color=color, lw=0.6,
                     label=f"{label} − RFM")
    axes[0].set_ylabel("Brightness temperature [K]")
    axes[0].set_title(f"{geometry.capitalize()} brightness temperature", loc="left")
    axes[0].legend(frameon=False, ncol=4, fontsize=8)
    axes[1].set_ylabel("JURASSIC − RFM [K]")
    axes[1].set_title(f"{geometry.capitalize()} brightness temperature difference", loc="left")
    axes[1].axhline(0, color="0.25", lw=0.5)
    axes[1].legend(frameon=False, ncol=3, fontsize=8)
    axes[1].set_xlabel("Wavenumber [cm$^{-1}$]")
    for ax in axes:
        ax.grid(alpha=0.18, lw=0.5)
    fig.savefig(output, dpi=220)
    plt.close(fig)


def write_metrics(root, output):
    rows = []
    for method, label, _ in METHODS:
        for ray, height in enumerate(HEIGHTS):
            data = limb_data(root, method)[ray]
            rows.append({"geometry": "limb", "case": f"{height} km geometric",
                         "method": label, "quantity": "absolute relative radiance difference",
                         "unit": "%", **limb_statistics(data)})
        for geometry in ("nadir", "zenith"):
            values = [row[3] for row in bt_data(root, geometry, method)]
            rows.append({"geometry": geometry, "case": geometry, "method": label,
                         "quantity": "absolute brightness temperature difference",
                         "unit": "K", **statistics(values),
                         "maximum_channel_cm-1": "", "reference_at_maximum": "",
                         "absolute_difference_at_maximum": ""})
    with output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys(), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return rows


def plot_limb_summary(metrics, output, plt):
    fig, ax = plt.subplots(figsize=(6.2, 4.4), constrained_layout=True)
    x_limb = list(range(4))
    for method, label, color in METHODS:
        rows = [row for row in metrics if row["geometry"] == "limb" and row["method"] == label]
        ax.plot(x_limb, [row["median"] for row in rows], marker="o", linestyle="-",
                color=color, label=f"{label} median")
        ax.plot(x_limb, [row["p95"] for row in rows], marker="^", linestyle="--",
                color=color, label=f"{label} 95th percentile")
    ax.set_xticks(x_limb, [f"{height} km" for height in HEIGHTS])
    ax.set_ylabel("Absolute relative difference [%]")
    ax.set_title("Limb radiance accuracy")
    ax.legend(frameon=False, fontsize=8, ncol=2)
    ax.grid(axis="y", alpha=0.2, lw=0.5)
    fig.savefig(output, dpi=220)
    plt.close(fig)


def plot_bt_summary(metrics, output, plt):
    fig, ax = plt.subplots(figsize=(6.2, 4.4), constrained_layout=True)
    x_bt = [0, 1]
    width = 0.25
    for offset, (method, label, color) in zip((-width, 0, width), METHODS):
        rows = [row for row in metrics if row["geometry"] in ("nadir", "zenith")
                and row["method"] == label]
        ax.bar([value + offset for value in x_bt], [row["rms"] for row in rows], width,
               color=color, label=f"{label} RMS")
        ax.scatter([value + offset for value in x_bt], [row["p95"] for row in rows],
                   color="black", marker="_", s=90, zorder=3)
    ax.set_xticks(x_bt, ("Nadir", "Zenith"))
    ax.set_ylabel("Absolute BT difference [K]")
    ax.set_title("Brightness temperature accuracy")
    ax.legend(frameon=False, fontsize=8, title="Black mark: 95th percentile")
    ax.grid(axis="y", alpha=0.2, lw=0.5)
    fig.savefig(output, dpi=220)
    plt.close(fig)


def plot_runtime_summary(root, output, plt):
    timings = {}
    for method in ("ega", "rfm"):
        directory = "test_ega" if method == "ega" else "rfm_reference"
        with (root / directory / "timings.csv").open(newline="") as stream:
            for row in csv.DictReader(stream):
                timings[row["geometry"], method] = row
    with (root / "test_cga" / "timings.csv").open(newline="") as stream:
        cga_timings = {row["geometry"]: row for row in csv.DictReader(stream)}
    geometries = ("limb", "nadir", "zenith")
    cases = ("Limb\n(4 rays)", "Nadir\n(1 ray)", "Zenith\n(1 ray)")
    x = list(range(len(cases)))
    ega_time = [float(timings[geometry, "ega"]["model_s"])
                for geometry in geometries]
    cga_time = [float(cga_timings[geometry]["model_s"])
                for geometry in geometries]
    rfm_time = [float(timings[geometry, "rfm"]["model_s"])
                for geometry in geometries]
    fig, ax = plt.subplots(figsize=(6.2, 4.4), constrained_layout=True)
    width = 0.25
    ax.bar([value - width for value in x], ega_time, width, label="EGA", color="#0072B2")
    ax.bar(x, cga_time, width, label="CGA", color="#E69F00")
    ax.bar([value + width for value in x], rfm_time, width, label="RFM", color="0.35")
    ax.set_yscale("log")
    ax.set_xticks(x, cases)
    ax.set_ylabel("Model time per validation case [s]")
    ax.set_title("Single-core model runtime")
    ax.legend(frameon=False)
    ax.grid(axis="y", alpha=0.2, lw=0.5)
    fig.savefig(output, dpi=220)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path,
                        help="output directory (default: projects/validation/analysis)")
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    output_dir = (args.output_dir or root / "analysis").expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        parser.error(f"matplotlib is required: {exc}")
    plt.rcParams.update({"font.size": 9, "axes.titlesize": 10, "savefig.facecolor": "white"})
    metrics = write_metrics(root, output_dir / "accuracy_metrics.csv")
    plot_limb_radiance(root, output_dir / "limb_radiance_spectra.png", plt)
    plot_limb_relative(root, output_dir / "limb_relative_errors.png", plt)
    plot_bt(root, "nadir", output_dir / "nadir_brightness_temperature.png", plt)
    plot_bt(root, "zenith", output_dir / "zenith_brightness_temperature.png", plt)
    plot_limb_summary(metrics, output_dir / "limb_accuracy_summary.png", plt)
    plot_bt_summary(metrics, output_dir / "brightness_temperature_accuracy_summary.png", plt)
    plot_runtime_summary(root, output_dir / "runtime_summary.png", plt)
    print(f"Metrics and figures: {output_dir}")


if __name__ == "__main__":
    main()
