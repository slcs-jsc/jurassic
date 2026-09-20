#!/usr/bin/env python3
"""Generate channel-matched JURASSIC and optional RFM validation spectra."""
import argparse
import csv
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import math
import os
import re
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
# Gas order and cross-section mapping from the full_spec_v6 control template.
GASES = (
    "CO2", "H2O", "N2", "O2", "O3", "C2H2", "C2H6", "CCl4", "CH4",
    "ClO", "ClONO2", "CO", "COF2", "F11", "F12", "F14", "F22", "H2O2",
    "HCN", "HNO3", "HNO4", "N2O", "N2O5", "NH3", "NO", "NO2", "OCS",
    "SF6", "SO2", "F113", "F114", "HF", "HCl", "HBr", "HOCl", "H2CO",
)
XSC_FILES = {
    7: "CCl4.xsc", 10: "ClONO2.xsc", 13: "CFC-11.xsc",
    14: "CFC-12.xsc", 15: "CFC-14.xsc", 16: "HCFC-22.xsc",
    20: "HNO4.xsc", 22: "N2O5.xsc", 27: "SF6.xsc",
    29: "CFC-113.xsc", 30: "CFC-114.xsc",
}
GEOMETRY = {
    "limb": ("OBSZ = 780",),
    "nadir": ("OBSZ = 700", "LAT0 = 0", "LAT1 = 0"),
    "zenith": ("OBSZ = 0", "VPZ = 700", "THETA0 = 0", "THETA1 = 0"),
}
# TIMER_* values come from JURASSIC; RFM_PHASE_* values require the separate
# instrumented RFM executable. Keep both timing boundaries visible in outputs.
LIMB_HEIGHTS = (5, 10, 20, 50)
TIMER_NAMES = ("READ_CTL", "READ_TBL", "READ_TASK", "READ_DIRLIST",
               "READ_OBSREF", "READ_ATM", "READ_OBS", "FORMOD",
               "WRITE_OBS", "FINALIZE")


def hardware_metadata():
    """Record CPU identity and the logical CPUs available to this run."""
    cpu_model = "unknown"
    physical_cores = set()
    current_physical = None
    current_core = None
    cpuinfo = Path("/proc/cpuinfo")
    if cpuinfo.is_file():
        for line in [*cpuinfo.read_text().splitlines(), ""]:
            if line.startswith("model name") and cpu_model == "unknown":
                cpu_model = line.split(":", 1)[1].strip()
            elif line.startswith("physical id"):
                current_physical = line.split(":", 1)[1].strip()
            elif line.startswith("core id"):
                current_core = line.split(":", 1)[1].strip()
            elif not line.strip():
                if current_physical is not None and current_core is not None:
                    physical_cores.add((current_physical, current_core))
                current_physical = current_core = None
    affinity = sorted(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else []
    return {
        "cpu_model": cpu_model,
        "physical_cores": len(physical_cores) or None,
        "logical_cpus": os.cpu_count(),
        "process_affinity_logical_cpus": affinity,
    }


def sha256(path):
    """Return the SHA-256 digest of a potentially large external input."""
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def run(command, cwd, log, env):
    """Run one model command, capture its output, and return wall time."""
    start = time.perf_counter()
    result = subprocess.run(command, cwd=cwd, env=env, text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    elapsed = time.perf_counter() - start
    log.write_text(result.stdout)
    if result.returncode:
        raise RuntimeError(f"command failed ({result.returncode}): {' '.join(command)}; see {log}")
    return elapsed


def parse_timers(log):
    """Collect JURASSIC timers and sum all instrumented RFM invocations."""
    text = log.read_text()
    values = {}
    for name in TIMER_NAMES:
        match = re.search(rf"^TIMER_{name} = .*?mean= ([0-9.eE+-]+) s", text, re.MULTILINE)
        if not match:
            raise ValueError(f"missing TIMER_{name} in {log}")
        values[name.lower() + "_s"] = float(match.group(1))
    values["input_setup_s"] = sum(values[name.lower() + "_s"] for name in TIMER_NAMES[:7])
    values["forward_section_s"] = values["formod_s"]
    values["output_s"] = values["write_obs_s"]
    # FORMOD=2 launches RFM once per overlapping spectral block. Sum every
    # invocation, not just the first timer line.
    for prefix, names in (("RFM_PHASE_", ("DRIVER_S", "PROFILE_S", "PATH_S", "SPECTRAL_S", "OUTPUT_S", "SPECTRAL_EX_OUTPUT_S")),
                          ("RFM_HITRAN_", ("BIN_READ_S", "BIN_READ_CALLS", "INIT_S", "INIT_CALLS"))):
        for name in names:
            matches = re.findall(rf"^{prefix}{name}=\s*([0-9.eE+-]+)$", text, re.MULTILINE)
            if matches:
                value = sum(float(item) for item in matches)
                key = ("rfm_phase_" if prefix == "RFM_PHASE_" else "rfm_hitran_") + name.lower()
                values[key] = int(value) if name.endswith("CALLS") else value
    return values


def read_science(path, nus):
    """Read channel radiances from one JURASSIC NetCDF result."""
    import netCDF4
    with netCDF4.Dataset(path) as dataset:
        rows = None
        for nu in nus:
            name = f"rad_{nu:.4f}"
            if name not in dataset.variables:
                raise ValueError(f"missing radiance channel {name} in {path}")
            data = dataset.variables[name][:]
            if data.shape[0] != 1:
                raise ValueError(f"expected one profile in {path}")
            if rows is None:
                rows = [[] for _ in range(data.shape[1])]
            elif data.shape[1] != len(rows):
                raise ValueError(f"ray count differs by channel in {path}")
            for ray, value in enumerate(data[0]):
                if not math.isfinite(value):
                    raise ValueError(f"non-finite radiance in {path}, {name}, ray {ray}")
                rows[ray].append(float(value))
    if not rows:
        raise ValueError(f"empty spectrum: {path}")
    return rows


def merge_limb_observations(paths, view_heights, output):
    """Combine four one-ray ASCII observations in geometric-height order."""
    header = []
    rows = []
    for path, height in zip(paths, view_heights):
        lines = path.read_text().splitlines()
        if not header:
            header = [line for line in lines if line.startswith("#")]
        data = [line for line in lines if line.strip() and not line.startswith("#")]
        if len(data) != 1 or abs(float(data[0].split()[4]) - height) > 1e-3:
            raise ValueError(f"invalid limb ray for {height} km: {path}")
        rows.append(data[0])
    output.write_text("\n".join([*header, "", *rows]) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--tbl-dir", type=Path, default=os.environ.get("JURASSIC_TBL_DIR"),
                        help="external directory containing tria_<gas>.nc for all 36 gases")
    parser.add_argument("--bin-dir", type=Path, default=ROOT / "src",
                        help="directory containing the JURASSIC executables")
    parser.add_argument("--rfm-bin", type=Path, help="external RFM executable; requires --rfm-hit")
    parser.add_argument("--rfm-hit", type=Path, help="external HITRAN line file used by RFM")
    parser.add_argument("--rfm-xsc-dir", type=Path, default=os.environ.get("RFM_XSC_DIR"),
                        help="external directory containing required RFM cross sections")
    parser.add_argument("--method", choices=("ega", "cga", "rfm"), default="ega",
                        help="single forward-model method to calculate")
    parser.add_argument("--require-rfm-timers", action="store_true",
                        help="fail if RFM source phase and HITRAN input timers are missing")
    parser.add_argument("--geometries", nargs="+", choices=GEOMETRY, default=list(GEOMETRY),
                        help="viewing geometries to calculate")
    parser.add_argument("--nu-start", type=int, default=500,
                        help="first channel center in cm^-1")
    parser.add_argument("--nu-end", type=int, default=2999,
                        help="last channel center in cm^-1, inclusive")
    parser.add_argument("--chunk-size", type=int, default=128,
                        help="channels per independent model process; maximum 128")
    parser.add_argument("--jobs", type=int, default=1, help="spectral chunks to run concurrently; each uses one CPU thread")
    parser.add_argument("--output", type=Path, default=ROOT / "projects/validation/work/common",
                        help="internal work/output directory")
    if len(sys.argv) == 1:
        parser.print_help()
        return
    args = parser.parse_args()
    if args.nu_start > args.nu_end or args.chunk_size < 1 or args.chunk_size > 128 or args.jobs < 1:
        parser.error("require nu-start <= nu-end, 1 <= chunk-size <= 128, and jobs >= 1")
    if args.method == "rfm" and os.environ.get("OMP_NUM_THREADS", "1") != "1":
        parser.error("RFM comparison requires OMP_NUM_THREADS=1 for a single-core comparison")
    if args.tbl_dir is None:
        parser.error("set --tbl-dir or JURASSIC_TBL_DIR")
    table_dir = args.tbl_dir.expanduser().resolve()
    if not table_dir.is_dir():
        parser.error(f"lookup-table directory does not exist: {table_dir}")
    if args.method == "rfm" and (not args.rfm_bin or not args.rfm_hit):
        parser.error("RFM method requires --rfm-bin and --rfm-hit")
    if args.method != "rfm" and (args.rfm_bin or args.rfm_hit):
        parser.error("RFM inputs are only valid with --method rfm")
    if args.require_rfm_timers and not args.rfm_bin:
        parser.error("--require-rfm-timers requires --rfm-bin")
    rfm_bin = args.rfm_bin.expanduser().resolve() if args.rfm_bin else None
    rfm_hit = args.rfm_hit.expanduser().resolve() if args.rfm_hit else None
    for path in (rfm_bin, rfm_hit):
        if path and not path.is_file():
            parser.error(f"required RFM file missing: {path}")
    if rfm_bin and args.rfm_xsc_dir is None:
        parser.error("set --rfm-xsc-dir or RFM_XSC_DIR for the 36-gas RFM run")
    xsc_dir = args.rfm_xsc_dir.expanduser().resolve() if args.rfm_xsc_dir else None
    if rfm_bin:
        if not xsc_dir.is_dir():
            parser.error(f"RFM cross-section directory does not exist: {xsc_dir}")
        missing_xsc = [xsc_dir / name for name in XSC_FILES.values()
                       if not (xsc_dir / name).is_file()]
        if missing_xsc:
            parser.error(f"missing RFM cross sections: {', '.join(map(str, missing_xsc))}")
    if rfm_bin and not os.access(rfm_bin, os.X_OK):
        parser.error(f"RFM binary is not executable: {rfm_bin}")
    bins = {name: args.bin_dir.expanduser().resolve() / name for name in ("climatology", "formod", "obsfmt", "raytrace", *GEOMETRY)}
    for name in ("climatology", "formod", "obsfmt", *(("raytrace",) if "limb" in args.geometries else ()), *args.geometries):
        if not bins[name].is_file():
            parser.error(f"missing JURASSIC executable: {bins[name]}")
    nus = list(range(args.nu_start, args.nu_end + 1))
    try:
        import netCDF4
    except ImportError:
        parser.error("Python netCDF4 is required to check TRIA channel coverage")
    for gas in GASES:
        path = table_dir / f"tria_{gas}.nc"
        if not path.is_file():
            parser.error(f"missing lookup table: {path}")
        if gas in ("CO2", "H2O"):
            with netCDF4.Dataset(path) as table:
                missing = [nu for nu in nus if f"tbl_{nu:.4f}" not in table.variables]
            if missing:
                parser.error(f"{path} lacks {len(missing)} requested channels; first missing: {missing[0]} cm^-1")
    # The manifest records external inputs and settings, but a Git commit alone
    # does not capture uncommitted changes to this runner or the executables.
    execution_start = time.perf_counter()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    env = {**os.environ, "LC_ALL": "C", "OMP_NUM_THREADS": os.environ.get("OMP_NUM_THREADS", "1")}
    methods = {args.method: {"ega": "1", "cga": "0", "rfm": "2"}[args.method]}
    manifest = {
        "git_commit": subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True,
                                     capture_output=True, check=True).stdout.strip(),
        "table_dir": str(table_dir), "table_files": {gas: str(table_dir / f"tria_{gas}.nc") for gas in GASES},
        "gases": GASES,
        "cross_sections": {GASES[index]: {"path": str(xsc_dir / name), "sha256": sha256(xsc_dir / name)}
                           for index, name in XSC_FILES.items()} if rfm_bin else None,
        "nu_start": args.nu_start, "nu_end": args.nu_end,
        "chunk_size": args.chunk_size, "parallel_jobs": args.jobs, "geometries": args.geometries,
        "limb_geometric_tangent_heights_km": LIMB_HEIGHTS,
        "methods": methods, "rfm_binary": str(rfm_bin) if rfm_bin else None,
        "rfm_hitran": str(rfm_hit) if rfm_hit else None,
        "rfm_hitran_sha256": sha256(rfm_hit) if rfm_hit else None,
        "rfm_binary_sha256": sha256(rfm_bin) if rfm_bin else None,
        "rfm_source_timers_required": args.require_rfm_timers,
        "omp_num_threads": env["OMP_NUM_THREADS"],
        "hardware": hardware_metadata(),
        "timing_definition": "sum of per-chunk formod wall and named timer sections; excludes input generation and plots; chunks may run concurrently",
        "model_time_definition": {"jurassic": "sum of TIMER_FORMOD across chunks",
                                  "rfm": "sum of RFM PATH + SPECTRAL - OUTPUT across spectral-block invocations and chunks",
                                  "per_spectrum": "geometry model time divided by number of rays; includes all requested channels",
                                  "rfm_minus_measured_hitran_input": "RFM model time minus timed HITRAN initialization and binary record READs"},
        "rfm_forward_section_includes_internal_io": True,
        "status": "running",
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    summary = []
    for geometry in args.geometries:
        geom_dir = output / geometry
        geom_dir.mkdir(exist_ok=True)
        values = {method: [] for method in methods}
        seconds = {method: 0.0 for method in methods}
        phase_totals = {method: {} for method in methods}
        rfm_phase_complete = True
        chunks = [nus[index:index + args.chunk_size]
                  for index in range(0, len(nus), args.chunk_size)]

        def process_chunk(channels):
            # Every worker has its own directory because FORMOD=2 writes RFM
            # driver, spectrum, and log files in its current working directory.
            chunk = geom_dir / f"{channels[0]:04d}_{channels[-1]:04d}"
            chunk.mkdir(exist_ok=True)
            ctl = chunk / "case.ctl"
            lines = ["# Generated validation control", f"TBLBASE = {table_dir / 'tria'}", "TBLFMT = 3",
                     "WRITE_BBT = 0", "OBSFMT = 3", "REFRAC = 1", f"NG = {len(GASES)}", *[f"EMITTER[{i}] = {gas}" for i, gas in enumerate(GASES)],
                     f"ND = {len(channels)}", *[f"NU[{i}] = {nu:.4f}" for i, nu in enumerate(channels)],
                     *GEOMETRY[geometry]]
            ctl.write_text("\n".join(lines) + "\n")
            atm, obs = chunk / "atm.tab", chunk / "obs.nc"
            climate_command = [str(bins["climatology"]), str(ctl), str(atm),
                               "Z0", "0", "Z1", "90", "DZ", "1", "ZSURF", "1"]
            run(climate_command, chunk, chunk / "climatology.log", env)
            if geometry == "limb":
                # Z0/Z1 set geometric tangent heights. Both forward models use
                # refraction; RFM receives JURASSIC's traced tangent heights.
                geometry_commands = []
                convert_commands = []
                trace_commands = []
                ray_files = []
                for height in LIMB_HEIGHTS:
                    ray_file = chunk / f"ray_{height}km.tab"
                    command = [str(bins["limb"]), str(ctl), str(ray_file),
                               "OBSFMT", "1", "Z0", str(height), "Z1", str(height)]
                    run(command, chunk, chunk / f"geometry_{height}km.log", env)
                    geometry_commands.append(command)
                    ray_files.append(ray_file)
                combined = chunk / "obs.tab"
                merge_limb_observations(ray_files, LIMB_HEIGHTS, combined)
                convert_command = [str(bins["obsfmt"]), str(ctl), str(combined), "1", str(obs), "3"]
                run(convert_command, chunk, chunk / "obsfmt.log", env)
                convert_commands.append(convert_command)
                trace = chunk / "raytrace.tab"
                trace_command = [str(bins["raytrace"]), str(ctl), str(obs), str(atm), str(trace)]
                run(trace_command, chunk, chunk / "raytrace.log", env)
                trace_commands.append(trace_command)
                measured = [float(line.split()[7]) for line in trace.read_text().splitlines()
                            if line.strip() and not line.startswith("#")]
                if len(measured) != len(LIMB_HEIGHTS):
                    raise ValueError(f"expected four traced limb rays in {chunk}")
                (chunk / "limb_heights.json").write_text(json.dumps(
                    {"geometric_tangent_heights_km": LIMB_HEIGHTS,
                     "refracted_tangent_heights_km": measured}, indent=2) + "\n")
            else:
                geometry_commands = [[str(bins[geometry]), str(ctl), str(obs)]]
                convert_commands = []
                trace_commands = []
                run(geometry_commands[0], chunk, chunk / "geometry.log", env)
            (chunk / "commands.json").write_text(json.dumps(
                {"climatology": climate_command, "geometry": geometry_commands,
                 "obsfmt": convert_commands, "raytrace": trace_commands}, indent=2) + "\n")
            result_by_method = {}
            for method, code in methods.items():
                method_dir = chunk / method
                method_dir.mkdir(exist_ok=True)
                result = method_dir / "rad.nc"
                command = [str(bins["formod"]), str(ctl), str(obs), str(atm), str(result), "FORMOD", code]
                if method == "rfm":
                    # FORMOD=2 calls RFM and applies the TRIA channel response
                    # before writing radiances, so comparisons are channel matched.
                    command.extend(["RFMBIN", str(rfm_bin), "RFMHIT", str(rfm_hit)])
                    for gas_index, name in XSC_FILES.items():
                        command.extend([f"RFMXSC[{gas_index}]", str(xsc_dir / name)])
                (method_dir / "command.json").write_text(json.dumps(command, indent=2) + "\n")
                log = method_dir / "formod.log"
                wall = run(command, method_dir, log, env)
                timing = parse_timers(log)
                timer_complete = True
                if method == "rfm":
                    required = ("rfm_phase_driver_s", "rfm_phase_path_s",
                                "rfm_phase_spectral_s", "rfm_phase_output_s",
                                "rfm_hitran_init_s", "rfm_hitran_bin_read_s")
                    missing = [name for name in required if name not in timing]
                    if missing:
                        timer_complete = False
                        if args.require_rfm_timers:
                            raise ValueError(f"RFM source timers missing in {log}: {', '.join(missing)}")
                timing["wall_s"] = wall
                (method_dir / "timing.json").write_text(json.dumps(timing, indent=2) + "\n")
                result_by_method[method] = (wall, timing, result, timer_complete)
                print(f"{geometry} {chunk.name} {method}: complete ({wall:.1f} s)", flush=True)
            return result_by_method

        # map returns chunks in spectral order, independent of completion order.
        # The per-method sums remain comparable single-thread process times;
        # parallel elapsed job time is a separate quantity.
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            for channels, result_by_method in zip(chunks, pool.map(process_chunk, chunks)):
                for method, (wall, timing, result, timer_complete) in result_by_method.items():
                    rows = read_science(result, channels)
                    seconds[method] += wall
                    if method == "rfm" and not timer_complete:
                        rfm_phase_complete = False
                    for key, value in timing.items():
                        phase_totals[method][key] = phase_totals[method].get(key, 0) + value
                    if values[method] and len(rows) != len(values[method]):
                        raise ValueError(f"ray count changed between chunks for {geometry}")
                    if not values[method]:
                        values[method] = [[] for _ in rows]
                    for ray, row in enumerate(rows):
                        values[method][ray].extend(row)
        for method, spectra in values.items():
            directory = geom_dir / method
            directory.mkdir(exist_ok=True)
            with (directory / "spectrum.csv").open("w", newline="") as stream:
                writer = csv.writer(stream, lineterminator="\n")
                writer.writerow(("nu_cm-1", "ray", "radiance_W_m-2_sr-1_cm"))
                writer.writerows((nu, ray, value) for ray, row in enumerate(spectra)
                                 for nu, value in zip(nus, row))
            phases = phase_totals[method]
            phases["other_wall_s"] = max(0.0, phases["wall_s"] -
                phases["input_setup_s"] - phases["forward_section_s"] -
                phases["output_s"] - phases["finalize_s"])
            phases["spectrum_count"] = len(spectra)
            phases["channel_count"] = len(nus)
            if method == "rfm":
                if rfm_phase_complete:
                    # Remove RFM output time from its PATH + SPECTRAL phases;
                    # this still includes HITRAN input and spectral setup.
                    phases["model_s"] = (phases["rfm_phase_path_s"] +
                                         phases["rfm_phase_spectral_s"] -
                                         phases["rfm_phase_output_s"])
            else:
                phases["model_s"] = phases["formod_s"]
            if "model_s" in phases:
                # Limb has multiple rays per call; this is an amortized time per
                # complete spectrum, not single-ray latency.
                phases["model_s_per_spectrum"] = phases["model_s"] / len(spectra)
            if method == "rfm" and "model_s" in phases and all(
                    name in phases for name in ("rfm_hitran_init_s", "rfm_hitran_bin_read_s")):
                phases["model_minus_timed_hitran_input_s"] = (
                    phases["model_s"] - phases["rfm_hitran_init_s"] -
                    phases["rfm_hitran_bin_read_s"])
                phases["model_minus_timed_hitran_input_s_per_spectrum"] = (
                    phases["model_minus_timed_hitran_input_s"] / len(spectra))
            summary.append({"geometry": geometry, "method": method, **phases})
        print(f"{geometry}: {args.method.upper()} spectrum generated")
    with (output / "timings.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=("geometry", "method", "wall_s",
            *(name.lower() + "_s" for name in TIMER_NAMES), "input_setup_s",
            "forward_section_s", "output_s", "other_wall_s",
            "spectrum_count", "channel_count", "model_s", "model_s_per_spectrum",
            "model_minus_timed_hitran_input_s",
            "model_minus_timed_hitran_input_s_per_spectrum",
            "rfm_phase_driver_s", "rfm_phase_profile_s", "rfm_phase_path_s",
            "rfm_phase_spectral_s", "rfm_phase_output_s",
            "rfm_phase_spectral_ex_output_s", "rfm_hitran_bin_read_s",
            "rfm_hitran_bin_read_calls", "rfm_hitran_init_s",
            "rfm_hitran_init_calls"), lineterminator="\n")
        writer.writeheader()
        writer.writerows(summary)
    manifest["status"] = "complete"
    manifest["execution_wall_s"] = time.perf_counter() - execution_start
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Results: {output}")


if __name__ == "__main__":
    main()
