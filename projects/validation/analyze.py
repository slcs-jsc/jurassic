#!/usr/bin/env python3
"""Create geometry-appropriate accuracy plots from a completed validation run."""

import argparse
import csv
import json
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
            values = [float(row["relative_difference_percent"])
                      for row in limb_data(root, method)[ray]]
            rows.append({"geometry": "limb", "case": f"{height} km geometric",
                         "method": label, "quantity": "absolute relative radiance difference",
                         "unit": "%", **statistics(values)})
        for geometry in ("nadir", "zenith"):
            values = [row[3] for row in bt_data(root, geometry, method)]
            rows.append({"geometry": geometry, "case": geometry, "method": label,
                         "quantity": "absolute brightness temperature difference",
                         "unit": "K", **statistics(values)})
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


def write_report(root, metrics, output):
    """Write a self-contained Markdown overview of the current results."""
    with (root / "rfm_reference" / "manifest.json").open() as stream:
        manifest = json.load(stream)

    hardware = manifest.get("hardware", {})
    affinity = hardware.get("process_affinity_logical_cpus", [])
    affinity_text = ", ".join(str(cpu) for cpu in affinity) if affinity else "not recorded"

    timings = {}
    for method, directory in (("EGA", "test_ega"), ("CGA", "test_cga"),
                              ("RFM", "rfm_reference")):
        with (root / directory / "timings.csv").open(newline="") as stream:
            for row in csv.DictReader(stream):
                timings[row["geometry"], method] = float(row["model_s"])

    lines = [
        "# JURASSIC–RFM validation report",
        "",
        "This report is generated by `../analyze.py` from the compact spectra and timing",
        "files in this validation project. No values below are entered manually.",
        "",
        "## Configuration",
        "",
        f"- JURASSIC commit: `{manifest['git_commit']}`",
        f"- Spectral grid: {manifest['nu_start']}–{manifest['nu_end']} cm⁻¹ at 1 cm⁻¹ sampling",
        f"- Atmospheric composition: mid-latitude climatology with {len(manifest['gases'])} gases",
        "- Geometries: limb at 5, 10, 20, and 50 km geometric tangent height; one nadir and one zenith ray",
        "- Refraction: enabled consistently for JURASSIC and RFM",
        "- JURASSIC modes: EGA and CGA",
        "- Threads per model process: 1",
        "- Channels compared per spectrum: 2500",
        "- Spectral execution: 20 chunks of at most 128 channels; one contiguous RFM block per chunk",
        "",
        "RFM spectra are averaged with the same channel response functions used by",
        "JURASSIC. Limb errors are relative radiance errors. Nadir and zenith errors",
        "are absolute brightness temperature errors. All limb channels are included;",
        "only an exactly zero RFM radiance would have an undefined relative error.",
        "",
        "## Interpretation of the approximation errors",
        "",
        "JURASSIC replaces monochromatic radiative transfer by channel-averaged",
        "emissivities and Planck functions. EGA follows emissivity growth along an",
        "inhomogeneous ray path, whereas CGA replaces that path by an equivalent",
        "homogeneous path. Both remain band approximations: correlations between the",
        "spectral variation of the Planck function and emissivity along the path are",
        "not represented exactly. In addition, JURASSIC combines the channel-mean",
        "transmissions of individual gases multiplicatively. This neglects spectral",
        "correlation terms between overlapping absorption structures of different",
        "gases. These mechanisms are described for JURASSIC by",
        "[Baumeister and Hoffmann (2022)](https://doi.org/10.5194/gmd-15-1855-2022).",
        "",
        "The largest errors below occur in individual channels and should be read",
        "together with the median and 95th-percentile statistics. Published accuracies",
        "are specific to their setup: Gordley and Russell (1981) reported about 0.5%",
        "for a single-gas, approximately 100 cm^-1 broadband limb calculation, while",
        "Francis et al. (2006) reported channel-dependent radiance accuracies of",
        "0.5-1.0% or better for the 21-channel HIRDLS fast model, which combined CGA,",
        "EGA, and statistical regression. Neither result is a universal bound for the",
        "present 1 cm^-1, 36-gas calculation. Errors depend on spectral interval,",
        "channel response, atmospheric state, gas overlap, and viewing geometry.",
        "Consequently, each instrument or application requires its own line-by-line",
        "validation with the applicable spectral response functions.",
        "",
        "References: [Gordley and Russell (1981)](https://doi.org/10.1364/AO.20.000807);",
        "[Marshall et al. (1994)](https://doi.org/10.1016/0022-4073(94)90026-4);",
        "[Francis et al. (2006)](https://doi.org/10.1029/2005JD006270).",
        "",
        "## Limb spectra and errors",
        "",
        "![Limb radiance spectra](limb_radiance_spectra.png)",
        "",
        "![Relative limb radiance differences](limb_relative_errors.png)",
        "",
        "![Limb accuracy summary](limb_accuracy_summary.png)",
        "",
        "### Limb error statistics",
        "",
        "| Height | Method | RMS [%] | Median [%] | 95th percentile [%] | Maximum [%] |",
        "|---:|:---|---:|---:|---:|---:|",
    ]
    for row in metrics:
        if row["geometry"] == "limb":
            lines.append(f"| {row['case']} | {row['method']} | {row['rms']:.3f} | "
                         f"{row['median']:.3f} | {row['p95']:.3f} | {row['maximum']:.3f} |")

    lines += [
        "",
        "## Nadir and zenith brightness temperature",
        "",
        "![Nadir brightness temperature](nadir_brightness_temperature.png)",
        "",
        "![Zenith brightness temperature](zenith_brightness_temperature.png)",
        "",
        "![Brightness temperature accuracy summary](brightness_temperature_accuracy_summary.png)",
        "",
        "### Brightness temperature error statistics",
        "",
        "| Geometry | Method | RMS [K] | Median [K] | 95th percentile [K] | Maximum [K] |",
        "|:---|:---|---:|---:|---:|---:|",
    ]
    for row in metrics:
        if row["geometry"] in ("nadir", "zenith"):
            lines.append(f"| {row['geometry'].capitalize()} | {row['method']} | "
                         f"{row['rms']:.3f} | {row['median']:.3f} | "
                         f"{row['p95']:.3f} | {row['maximum']:.3f} |")

    lines += [
        "",
        "## Runtime",
        "",
        "![Single-core model runtime](runtime_summary.png)",
        "",
        "Reference timing hardware and execution:",
        "",
        f"- Processor: {hardware.get('cpu_model', 'not recorded')}",
        f"- CPU topology: {hardware.get('physical_cores', 'unknown')} physical cores, "
        f"{hardware.get('logical_cpus', 'unknown')} logical CPUs",
        f"- Execution: {manifest.get('parallel_jobs', 'unknown')} concurrent single-thread "
        f"processes restricted to logical CPUs {affinity_text}",
        "",
        "The model times exclude validation input generation, plotting, and final output",
        "writing. JURASSIC time is `TIMER_FORMOD`. RFM time is its measured path plus",
        "spectral phases minus measured output time. Times are totals for each validation",
        "case: limb contains four jointly calculated rays, while nadir and zenith contain",
        "one ray each. RFM shares spectral setup and HITRAN processing across the four",
        "limb rays; this is one joint calculation rather than four independent runs.",
        "Lookup-table reading and preparation are not included.",
        "",
        "| Geometry | EGA [s] | CGA [s] | RFM [s] | RFM/EGA | RFM/CGA |",
        "|:---|---:|---:|---:|---:|---:|",
    ]
    for geometry in ("limb", "nadir", "zenith"):
        ega = timings[geometry, "EGA"]
        cga = timings[geometry, "CGA"]
        rfm = timings[geometry, "RFM"]
        label = "Limb (4 rays)" if geometry == "limb" else f"{geometry.capitalize()} (1 ray)"
        lines.append(f"| {label} | {ega:.2f} | {cga:.2f} | {rfm:.2f} | "
                     f"{rfm / ega:.0f}× | {rfm / cga:.0f}× |")

    lines += [
        "",
        "These timings describe this recorded run and machine; they are not portable",
        "performance guarantees. Accuracy statistics are computed from the complete",
        "500–2999 cm⁻¹ channel set stored in the repository.",
        "",
    ]
    output.write_text("\n".join(lines))


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
    write_report(root, metrics, output_dir / "REPORT.md")
    print(f"Figures and report: {output_dir}")


if __name__ == "__main__":
    main()
