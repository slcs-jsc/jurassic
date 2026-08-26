#!/usr/bin/env python3
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

from validationlib import (
    REPO_ROOT,
    ROOT,
    render_ctl,
    select_cases,
    sha256_file,
    write_json,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate immutable JURASSIC validation references.")
    parser.add_argument("--profile", choices=("smoke", "full"))
    parser.add_argument("--case", action="append", default=[], dest="cases")
    parser.add_argument("--tblbase", default=os.environ.get("VALIDATION_TBLBASE"))
    parser.add_argument("--references", type=Path, default=ROOT / "references")
    parser.add_argument("--replace", action="store_true", help="Replace explicitly selected existing references.")
    parser.add_argument("--missing", action="store_true", help="Skip references that already exist.")
    parser.add_argument("--allow-dirty", action="store_true")
    return parser.parse_args()


def git_output(*args: str) -> str:
    return subprocess.check_output(("git", *args), cwd=REPO_ROOT, text=True).strip()


def run(command: list[str], log: Path | None = None) -> None:
    if log:
        with log.open("w") as handle:
            subprocess.run(command, cwd=log.parent, stdout=handle, stderr=subprocess.STDOUT, check=True)
    else:
        subprocess.run(command, check=True)


def main() -> int:
    args = parse_args()
    if not args.tblbase:
        print("VALIDATION_TBLBASE or --tblbase is required", file=sys.stderr)
        return 2
    try:
        cases = select_cases(args.profile, args.cases)
    except ValueError as exc:
        print(str(exc), file=sys.stderr)
        return 2

    dirty = bool(git_output("status", "--porcelain"))
    if dirty and not args.allow_dirty:
        print("refusing to generate references from a dirty worktree; use --allow-dirty explicitly", file=sys.stderr)
        return 2
    if args.missing and args.replace:
        print("--missing and --replace are mutually exclusive", file=sys.stderr)
        return 2
    commit = git_output("rev-parse", "HEAD")
    args.references.mkdir(parents=True, exist_ok=True)
    src = REPO_ROOT / "src"
    for executable in ("climatology", "formod"):
        if not (src / executable).is_file():
            print(f"missing executable: {src / executable}", file=sys.stderr)
            return 2

    for case in cases:
        geometry_exe = src / case.geometry
        if not geometry_exe.is_file():
            print(f"missing executable: {geometry_exe}", file=sys.stderr)
            return 2
        target = args.references / case.case_name
        if target.exists():
            if args.missing:
                print(f"skipped existing {target}")
                continue
            if not args.replace:
                print(f"reference already exists: {target}; use --replace explicitly", file=sys.stderr)
                return 2
            if not args.cases:
                print("--replace requires at least one explicit --case", file=sys.stderr)
                return 2

        with tempfile.TemporaryDirectory(prefix=f"{case.case_name}.", dir=args.references) as tmp_name:
            work = Path(tmp_name)
            canonical_ctl = work / "case.ctl"
            runtime_ctl = work / "runtime.ctl"
            canonical_ctl.write_text(render_ctl(case))
            runtime_ctl.write_text(render_ctl(case, str(Path(args.tblbase).expanduser().resolve())))
            run([str(src / "climatology"), str(runtime_ctl), "atm.tab"], work / "climatology.log")
            run([str(geometry_exe), str(runtime_ctl), "obs.nc"], work / "geometry.log")
            run(
                [str(src / "formod"), str(runtime_ctl), "obs.nc", "atm.tab", "rad.nc", "EXECUTION", "scalar"],
                work / "formod.log",
            )
            metadata = {
                "case": case.__dict__,
                "channels": len(case.channels),
                "created_utc": datetime.now(timezone.utc).isoformat(),
                "reference_git_commit": commit,
                "reference_worktree_dirty": dirty,
                "tblbase_at_generation": str(Path(args.tblbase).expanduser().resolve()),
                "files": {
                    name: sha256_file(work / name)
                    for name in ("case.ctl", "atm.tab", "obs.nc", "rad.nc")
                },
            }
            write_json(work / "metadata.json", metadata)
            runtime_ctl.unlink()
            if target.exists():
                shutil.rmtree(target)
            shutil.move(str(work), target)
        print(f"generated {target}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
