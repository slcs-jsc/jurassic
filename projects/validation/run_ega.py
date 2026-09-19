#!/usr/bin/env python3
"""Run the complete single-core JURASSIC EGA validation spectrum calculation."""

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
GEOMETRIES = ("limb", "nadir", "zenith")


def sha256(path):
    """Return the SHA-256 digest recorded in a result manifest."""
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def make_parser(description):
    """Create a parser that documents the fixed scientific setup."""
    return argparse.ArgumentParser(
        description=description,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog=("Fixed validation setup: 36 gases, 500--2999 cm^-1 at 1 cm^-1, "
                "limb 5/10/20/50 km plus nadir and zenith, TBLFMT 3, refraction enabled. "
                "Environment overrides: JURASSIC_TBL_DIR, VALIDATION_JOBS, "
                "VALIDATION_CPUSET."))


def add_common_arguments(parser):
    """Add execution options shared by the EGA, CGA, and RFM runners."""
    parser.add_argument("--tbl-dir", type=Path,
                        default=os.environ.get("JURASSIC_TBL_DIR", ROOT / "tab/tria_1cm/nc_1e-6"),
                        help="directory containing the external tria_<gas>.nc lookup tables")
    parser.add_argument("--bin-dir", type=Path, default=ROOT / "src",
                        help="directory containing the JURASSIC executables")
    parser.add_argument("--jobs", type=int, default=int(os.environ.get("VALIDATION_JOBS", "2")),
                        help="independent spectral chunks calculated concurrently")
    parser.add_argument("--cpuset", default=os.environ.get("VALIDATION_CPUSET", "0,2"),
                        help="taskset CPU list; use an empty value to disable CPU affinity")
    parser.add_argument("--force", action="store_true",
                        help="replace an existing completed result after the new run succeeds")


def write_spectra(work, method, output):
    """Merge per-geometry work files into one compact spectrum table."""
    with output.open("w", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(("geometry", "nu_cm-1", "ray", "radiance_W_m-2_sr-1_cm"))
        for geometry in GEOMETRIES:
            source = work / geometry / method / "spectrum.csv"
            with source.open(newline="") as input_stream:
                for row in csv.DictReader(input_stream):
                    writer.writerow((geometry, row["nu_cm-1"], row["ray"],
                                     row["radiance_W_m-2_sr-1_cm"]))


def geometry_rows(work):
    """Extract the common ray definitions from the first spectral chunk."""
    import netCDF4
    rows = []
    names = ("time", "obs_z", "obs_lon", "obs_lat", "vp_z", "vp_lon", "vp_lat")
    for geometry in GEOMETRIES:
        chunk = sorted((work / geometry).glob("[0-9]*_[0-9]*"))[0]
        with netCDF4.Dataset(chunk / "obs.nc") as dataset:
            for ray in range(dataset.dimensions["ray"].size):
                rows.append((geometry, ray, *(float(dataset[name][0, ray]) for name in names)))
    return rows


def write_inputs(work, output):
    """Record the atmosphere and geometry used to create the RFM reference."""
    output.mkdir()
    first_chunk = sorted((work / "limb").glob("[0-9]*_[0-9]*"))[0]
    shutil.copy2(first_chunk / "atm.tab", output / "atmosphere.tab")
    with (output / "geometry.csv").open("w", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(("geometry", "ray", "time_s", "observer_z_km", "observer_lon_deg",
                         "observer_lat_deg", "view_z_km", "view_lon_deg", "view_lat_deg"))
        writer.writerows(geometry_rows(work))


def check_reference_inputs(work):
    """Require JURASSIC runs to reproduce the recorded RFM input state."""
    reference = HERE / "rfm_reference/input"
    first_chunk = sorted((work / "limb").glob("[0-9]*_[0-9]*"))[0]
    if (first_chunk / "atm.tab").read_bytes() != (reference / "atmosphere.tab").read_bytes():
        raise RuntimeError("generated atmosphere differs from rfm_reference/input/atmosphere.tab")
    with (reference / "geometry.csv").open(newline="") as stream:
        expected = list(csv.reader(stream))[1:]
    actual = [[str(value) for value in row] for row in geometry_rows(work)]
    if len(expected) != len(actual):
        raise RuntimeError("generated geometry differs from RFM reference")
    for old, new in zip(expected, actual):
        if old[:2] != new[:2] or any(abs(float(a) - float(b)) > 1e-9
                                     for a, b in zip(old[2:], new[2:])):
            raise RuntimeError("generated geometry differs from RFM reference")


def execute(method, target_name, args, extra=()):
    """Run one method and atomically replace its compact result directory."""
    target = HERE / target_name
    if target.exists() and any(target.iterdir()) and not args.force:
        raise RuntimeError(f"{target} already contains results; pass --force to replace them")
    work_root = HERE / "work"
    work = work_root / method
    if work.exists():
        if not args.force:
            raise RuntimeError(f"unfinished work directory exists: {work}; inspect it or pass --force")
        shutil.rmtree(work)
    work_root.mkdir(exist_ok=True)
    command = [sys.executable, str(HERE / "common.py"), "--method", method,
               "--tbl-dir", str(args.tbl_dir), "--bin-dir", str(args.bin_dir),
               "--jobs", str(args.jobs), "--output", str(work), *extra]
    if args.cpuset:
        command = ["taskset", "-c", args.cpuset, *command]
    env = {**os.environ, "OMP_NUM_THREADS": "1"}
    subprocess.run(command, check=True, env=env)

    # Build a complete replacement beside the current result.  The existing
    # compact result remains usable if calculation or validation fails.
    staging = HERE / f".{target_name}.new"
    if staging.exists():
        shutil.rmtree(staging)
    staging.mkdir()
    write_spectra(work, method, staging / "spectra.csv")
    shutil.copy2(work / "timings.csv", staging / "timings.csv")
    manifest = json.loads((work / "manifest.json").read_text())
    manifest["method"] = method
    manifest["generation_command"] = command
    manifest["spectrum_sha256"] = sha256(staging / "spectra.csv")
    manifest.pop("methods", None)
    if method == "rfm":
        write_inputs(work, staging / "input")
        manifest["input_sha256"] = {
            path.name: sha256(path) for path in sorted((staging / "input").iterdir())
        }
    else:
        check_reference_inputs(work)
        manifest["reference_spectrum_sha256"] = sha256(HERE / "rfm_reference/spectra.csv")
    (staging / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    if target.exists():
        shutil.rmtree(target)
    staging.rename(target)
    shutil.rmtree(work)
    print(f"Results: {target}")


def main():
    parser = make_parser(__doc__)
    add_common_arguments(parser)
    if len(sys.argv) == 1:
        parser.print_help()
        return
    args = parser.parse_args()
    try:
        execute("ega", "test_ega", args)
    except RuntimeError as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    main()
