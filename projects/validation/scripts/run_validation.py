#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import subprocess
import sys

from validationlib import (
    ROOT,
    compare_netcdf,
    render_ctl,
    select_cases,
    sha256_file,
    write_json,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Validate a JURASSIC build against frozen references.")
    parser.add_argument("--profile", choices=("smoke", "full"), default="smoke")
    parser.add_argument("--case", action="append", default=[], dest="cases")
    parser.add_argument("--tblbase", default=os.environ.get("VALIDATION_TBLBASE"))
    parser.add_argument("--references", type=Path, default=ROOT / "references")
    parser.add_argument("--formod", type=Path, default=ROOT.parents[1] / "src" / "formod")
    parser.add_argument("--execution", choices=("batch", "scalar"), default="batch")
    parser.add_argument("--formod-method", choices=("cga", "ega", "rfm"), default="ega")
    parser.add_argument("--max-abs-rel-percent", type=float)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if not args.tblbase:
        print("VALIDATION_TBLBASE or --tblbase is required", file=sys.stderr)
        return 2
    if not args.formod.is_file():
        print(f"missing formod executable: {args.formod}", file=sys.stderr)
        return 2
    if args.max_abs_rel_percent is not None and (
        not math.isfinite(args.max_abs_rel_percent) or args.max_abs_rel_percent < 0
    ):
        print("--max-abs-rel-percent must be finite and non-negative", file=sys.stderr)
        return 2
    try:
        cases = select_cases(args.profile if not args.cases else None, args.cases)
    except ValueError as exc:
        print(str(exc), file=sys.stderr)
        return 2
    output = args.output or ROOT / "runs" / datetime.now().strftime("validation_%Y%m%d_%H%M%S")
    output.mkdir(parents=True, exist_ok=False)
    rows: list[dict[str, object]] = []
    failures = 0
    for case in cases:
        reference = args.references / case.case_name
        metadata_path = reference / "metadata.json"
        case_dir = output / case.case_name
        case_dir.mkdir()
        errors: list[str] = []
        if not metadata_path.is_file():
            errors.append(f"missing reference metadata: {metadata_path}")
        else:
            try:
                metadata = json.loads(metadata_path.read_text())
            except (json.JSONDecodeError, OSError) as exc:
                errors.append(f"invalid reference metadata: {exc}")
                metadata = {}
            if metadata.get("case") != case.__dict__:
                errors.append("reference metadata does not match manifest")
            frozen_files = metadata.get("files", {})
            for name in ("case.ctl", "atm.tab", "obs.nc", "rad.nc"):
                path = reference / name
                if not path.is_file():
                    errors.append(f"missing frozen file: {name}")
                elif frozen_files.get(name) != sha256_file(path):
                    errors.append(f"checksum mismatch for frozen file: {name}")
            canonical_ctl = reference / "case.ctl"
            if canonical_ctl.is_file() and canonical_ctl.read_text() != render_ctl(case):
                errors.append("reference control file does not match manifest")
        result = {"status": "FAIL", "channels": 0, "max_abs_rel_percent": 0.0, "errors": errors}
        if not errors:
            runtime_ctl = case_dir / "case.ctl"
            runtime_ctl.write_text(render_ctl(case, str(Path(args.tblbase).expanduser().resolve())))
            log = case_dir / "formod.log"
            rad = case_dir / "rad.nc"
            command = [
                str(args.formod.resolve()), str(runtime_ctl),
                str(reference / "obs.nc"), str(reference / "atm.tab"), str(rad),
                "EXECUTION", args.execution, "FORMOD", str({"cga": 0, "ega": 1, "rfm": 2}[args.formod_method]),
                "OBSREF", str(reference / "rad.nc"),
            ]
            completed = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            log.write_text(completed.stdout)
            if completed.returncode:
                errors.append(f"formod exited with status {completed.returncode}")
            threshold = args.max_abs_rel_percent
            if threshold is None:
                threshold = case.max_abs_rel_percent
            result = compare_netcdf(
                reference / "rad.nc", rad, len(case.channels), threshold
            )
            result["errors"] = errors + list(result["errors"])
            if result["errors"]:
                result["status"] = "FAIL"
        write_json(case_dir / "accuracy.json", result)
        failures += result["status"] != "PASS"
        rows.append({
            "case_name": case.case_name,
            "profile": case.profile,
            "geometry": case.geometry,
            "channels": len(case.channels),
            "execution": args.execution,
            "formod_method": args.formod_method,
            "status": result["status"],
            "max_abs_rel_percent": result["max_abs_rel_percent"],
            "errors": " | ".join(result["errors"]),
        })
        print(f"{result['status']}\t{case.case_name}\tmax_abs_rel={result['max_abs_rel_percent']}%")
    with (output / "summary.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    write_json(output / "run.json", {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "formod": str(args.formod.resolve()),
        "execution": args.execution,
        "formod_method": args.formod_method,
        "tblbase": str(Path(args.tblbase).expanduser().resolve()),
    })
    print(f"results: {output}")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
