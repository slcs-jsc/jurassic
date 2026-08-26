#!/usr/bin/env python3
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import csv
import hashlib
import json
import math
import re


ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = ROOT.parents[1]
CONFIGS = ROOT / "configs"
MANIFEST = CONFIGS / "validation_cases.tsv"
TBLBASE_TOKEN = "__VALIDATION_TBLBASE__"


@dataclass(frozen=True)
class ValidationCase:
    case_name: str
    profile: str
    geometry: str
    nu_min: int
    nu_max: int
    nu_step: int
    gas_set: str
    ray_set: str
    max_abs_rel_percent: float

    @property
    def channels(self) -> tuple[int, ...]:
        return tuple(range(self.nu_min, self.nu_max + 1, self.nu_step))


def load_cases() -> list[ValidationCase]:
    with MANIFEST.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    cases: list[ValidationCase] = []
    names: set[str] = set()
    for lineno, row in enumerate(rows, 2):
        case = ValidationCase(
            case_name=row["case_name"].strip(),
            profile=row["profile"].strip(),
            geometry=row["geometry"].strip(),
            nu_min=int(row["nu_min"]),
            nu_max=int(row["nu_max"]),
            nu_step=int(row["nu_step"]),
            gas_set=row["gas_set"].strip(),
            ray_set=row["ray_set"].strip(),
            max_abs_rel_percent=float(row["max_abs_rel_percent"]),
        )
        if not case.case_name or case.case_name in names:
            raise ValueError(f"invalid or duplicate case_name at {MANIFEST}:{lineno}")
        names.add(case.case_name)
        if case.profile not in {"smoke", "full"}:
            raise ValueError(f"invalid profile for {case.case_name}: {case.profile}")
        if case.geometry not in {"zenith", "nadir", "limb"}:
            raise ValueError(f"invalid geometry for {case.case_name}: {case.geometry}")
        if case.ray_set != "medium":
            raise ValueError(f"unsupported ray set for {case.case_name}: {case.ray_set}")
        if case.nu_step <= 0 or case.nu_min > case.nu_max:
            raise ValueError(f"invalid spectral range for {case.case_name}")
        if case.channels[-1] != case.nu_max:
            raise ValueError(f"spectral range does not include nu_max for {case.case_name}")
        if not math.isfinite(case.max_abs_rel_percent) or case.max_abs_rel_percent < 0:
            raise ValueError(f"invalid tolerance for {case.case_name}")
        gas_file = CONFIGS / "gas_sets" / f"{case.gas_set}.txt"
        if not gas_file.is_file():
            raise ValueError(f"missing gas set for {case.case_name}: {gas_file}")
        cases.append(case)
    return cases


def select_cases(profile: str | None, names: list[str]) -> list[ValidationCase]:
    cases = load_cases()
    known = {case.case_name for case in cases}
    unknown = sorted(set(names) - known)
    if unknown:
        raise ValueError(f"unknown validation case(s): {', '.join(unknown)}")
    if names:
        selected = [case for case in cases if case.case_name in names]
    elif profile:
        selected = [case for case in cases if case.profile == profile]
    else:
        selected = cases
    if not selected:
        raise ValueError("no validation cases selected")
    return selected


def read_gases(case: ValidationCase) -> list[str]:
    path = CONFIGS / "gas_sets" / f"{case.gas_set}.txt"
    return [
        line.strip()
        for line in path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]


def render_ctl(case: ValidationCase, tblbase: str = TBLBASE_TOKEN) -> str:
    gases = read_gases(case)
    lines = [
        "# Generated from projects/validation/configs/validation_cases.tsv.",
        f"TBLBASE = {tblbase}",
        "TBLFMT = 3",
        "OBSFMT = 3",
        "",
        f"NG = {len(gases)}",
    ]
    lines.extend(f"EMITTER[{idx}] = {gas}" for idx, gas in enumerate(gases))
    lines.extend(("", f"ND = {len(case.channels)}"))
    lines.extend(f"NU[{idx}] = {nu:.4f}" for idx, nu in enumerate(case.channels))
    lines.append("")
    if case.geometry == "zenith":
        lines.extend(("VPZ = 700", "OBSZ = 0", "THETA0 = -90.0", "THETA1 = 90.0", "DTHETA = 2.857142857142857"))
    elif case.geometry == "nadir":
        # Nadir: 64 rays over ~2200 km (19.8 degrees) centered at equator
        lines.extend(("OBSZ = 700", "LAT0 = -9.9", "LAT1 = 9.9", "DLAT = 0.314285714285714"))
    else:
        lines.extend(("OBSZ = 780", "Z0 = 5", "Z1 = 68", "DZ = 1"))

    # WRITE_BBT: Brightness temperature for zenith/nadir, radiance for limb
    write_bbt = 0 if case.geometry == "limb" else 1
    lines.extend(("", f"WRITE_BBT = {write_bbt}", ""))
    return "\n".join(lines)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()



EVAL_RE = re.compile(
    r"EVAL: nu=\s*(?P<nu>\S+) cm\^-1"
    r" \| MeanRE=\s*(?P<mre>\S+) %"
    r" \| SDRE=\s*(?P<sdre>\S+) %"
    r" \| MinRE=\s*(?P<minre>\S+) %"
    r" \| MaxRE=\s*(?P<maxre>\S+) %"
)


def evaluate_log(text: str, expected_channels: int, threshold: float) -> dict[str, object]:
    matches = list(EVAL_RE.finditer(text))
    errors: list[str] = []
    maximum = 0.0
    maximum_nu: float | None = None
    if len(matches) != expected_channels:
        errors.append(f"expected {expected_channels} EVAL records, found {len(matches)}")
    for match in matches:
        try:
            values = {key: float(match.group(key)) for key in ("nu", "mre", "sdre", "minre", "maxre")}
        except ValueError:
            errors.append(f"invalid numeric EVAL record: {match.group(0)}")
            continue
        if not all(math.isfinite(value) for value in values.values()):
            errors.append(f"non-finite EVAL record at nu={values['nu']}")
            continue
        channel_maximum = max(abs(values["minre"]), abs(values["maxre"]))
        if maximum_nu is None or channel_maximum > maximum:
            maximum = channel_maximum
            maximum_nu = values["nu"]
    if maximum > threshold:
        errors.append(f"maximum absolute relative error {maximum:.9g}% exceeds {threshold:.9g}%")
    return {
        "status": "FAIL" if errors else "PASS",
        "channels": len(matches),
        "max_abs_rel_percent": maximum,
        "max_error_nu_cm-1": maximum_nu,
        "threshold_percent": threshold,
        "errors": errors,
    }


def compare_netcdf(
    reference_path: Path,
    candidate_path: Path,
    expected_channels: int,
    max_abs_rel_percent: float,
    bt_abs_tolerance: float = 1e-10,
    tau_abs_tolerance: float = 1e-12,
) -> dict[str, object]:
    try:
        import netCDF4
        import numpy as np
    except ImportError as exc:
        raise RuntimeError("Python netCDF4 and numpy are required for validation") from exc

    errors: list[str] = []
    maximum_relative = 0.0
    maximum_absolute = 0.0
    maximum_variable: str | None = None
    compared_values = 0
    relative_tolerance = max_abs_rel_percent / 100.0
    with netCDF4.Dataset(reference_path) as reference, netCDF4.Dataset(candidate_path) as candidate:
        reference_names = sorted(
            name for name in reference.variables if name.startswith(("bt_", "rad_", "tau_"))
        )
        candidate_names = sorted(
            name for name in candidate.variables if name.startswith(("bt_", "rad_", "tau_"))
        )
        reference_radiances = [name for name in reference_names if name.startswith(("bt_", "rad_"))]
        if len(reference_radiances) != expected_channels:
            errors.append(
                f"expected {expected_channels} radiance variables, found {len(reference_radiances)}"
            )
        if reference_names != candidate_names:
            errors.append("reference and candidate science-variable sets differ")
        for name in sorted(set(reference_names) & set(candidate_names)):
            ref = np.ma.filled(reference.variables[name][:], np.nan).astype(float)
            test = np.ma.filled(candidate.variables[name][:], np.nan).astype(float)
            if ref.shape != test.shape:
                errors.append(f"shape mismatch for {name}: {ref.shape} != {test.shape}")
                continue
            ref_finite = np.isfinite(ref)
            test_finite = np.isfinite(test)
            if not np.array_equal(ref_finite, test_finite):
                errors.append(f"finite/NaN mask mismatch for {name}")
                continue
            finite = ref_finite & test_finite
            if not np.any(finite):
                continue
            ref_values = ref[finite]
            test_values = test[finite]
            diff = np.abs(test_values - ref_values)
            absolute_tolerance = bt_abs_tolerance if name.startswith(("bt_", "rad_")) else tau_abs_tolerance
            allowed = absolute_tolerance + relative_tolerance * np.abs(ref_values)
            if np.any(diff > allowed):
                worst = int(np.argmax(diff / np.maximum(allowed, np.finfo(float).tiny)))
                errors.append(
                    f"tolerance exceeded for {name}: abs={diff[worst]:.9g}, "
                    f"allowed={allowed[worst]:.9g}, ref={ref_values[worst]:.9g}"
                )
            local_absolute = float(np.max(diff))
            nonzero = ref_values != 0
            local_relative = (
                float(np.max(100.0 * diff[nonzero] / np.abs(ref_values[nonzero])))
                if np.any(nonzero) else 0.0
            )
            if local_absolute > maximum_absolute:
                maximum_absolute = local_absolute
            if maximum_variable is None or local_relative > maximum_relative:
                maximum_relative = local_relative
                maximum_variable = name
            compared_values += int(np.count_nonzero(finite))
    return {
        "status": "FAIL" if errors else "PASS",
        "channels": expected_channels,
        "compared_values": compared_values,
        "max_abs_rel_percent": maximum_relative,
        "max_absolute_difference": maximum_absolute,
        "max_error_variable": maximum_variable,
        "threshold_percent": max_abs_rel_percent,
        "bt_abs_tolerance": bt_abs_tolerance,
        "tau_abs_tolerance": tau_abs_tolerance,
        "errors": errors,
    }

def write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
