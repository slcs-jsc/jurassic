#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

from validationlib import ROOT, render_ctl, select_cases, sha256_file


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Audit validation-reference completeness and provenance.")
    parser.add_argument("--profile", choices=("smoke", "full"))
    parser.add_argument("--case", action="append", default=[], dest="cases")
    parser.add_argument("--references", type=Path, default=ROOT / "references")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        cases = select_cases(args.profile, args.cases)
    except ValueError as exc:
        print(str(exc), file=sys.stderr)
        return 2
    failures = 0
    for case in cases:
        directory = args.references / case.case_name
        errors: list[str] = []
        metadata_path = directory / "metadata.json"
        if not metadata_path.is_file():
            errors.append("missing metadata.json")
        else:
            try:
                metadata = json.loads(metadata_path.read_text())
            except (json.JSONDecodeError, OSError) as exc:
                errors.append(f"invalid metadata: {exc}")
                metadata = {}
            if metadata.get("case") != case.__dict__:
                errors.append("manifest row differs from frozen metadata")
            if metadata.get("channels") != len(case.channels):
                errors.append("channel count differs from manifest")
            files = metadata.get("files", {})
            for name in ("case.ctl", "atm.tab", "obs.nc", "rad.nc"):
                path = directory / name
                if not path.is_file():
                    errors.append(f"missing {name}")
                elif files.get(name) != sha256_file(path):
                    errors.append(f"checksum mismatch for {name}")
            canonical = directory / "case.ctl"
            if canonical.is_file() and canonical.read_text() != render_ctl(case):
                errors.append("canonical control file is stale")
        status = "FAIL" if errors else "PASS"
        print(f"{status}\t{case.case_name}" + (f"\t{' | '.join(errors)}" if errors else ""))
        failures += bool(errors)
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
