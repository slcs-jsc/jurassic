#!/usr/bin/env python3
"""Measure the current CPU forward model using its built-in TASK=t timer."""
import argparse
import csv
import os
from pathlib import Path
import re
import subprocess

ROOT = Path(__file__).resolve().parents[2]
TIMING = re.compile(r"RUNTIME: mean= ([\d.eE+-]+) s \| stddev= ([\d.eE+-]+) s \| min= ([\d.eE+-]+) s \| max= ([\d.eE+-]+) s")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bin", type=Path, default=ROOT / "src/formod")
    parser.add_argument("--tbl-dir", type=Path, default=os.environ.get("JURASSIC_TBL_DIR"),
                        help="directory containing tria_<emitter>.nc (or set JURASSIC_TBL_DIR)")
    parser.add_argument("--threads", type=int, nargs="+", default=[1, 2, 4])
    parser.add_argument("--cases", nargs="+", choices=["limb", "nadir", "zenith"], default=["limb", "nadir", "zenith"])
    parser.add_argument("--output", type=Path, default=Path("projects/benchmark/runs/cpu.tsv"))
    args = parser.parse_args()
    if any(n < 1 for n in args.threads):
        parser.error("thread counts must be positive")
    if args.tbl_dir is None:
        parser.error("set --tbl-dir or JURASSIC_TBL_DIR to the external lookup-table directory")
    tbl_dir = args.tbl_dir.expanduser().resolve()
    if not tbl_dir.is_dir():
        parser.error(f"lookup-table directory does not exist: {tbl_dir}")
    for case in args.cases:
        control = (ROOT / "projects/examples" / case / f"{case}.ctl").read_text()
        gases = [line.split("=", 1)[1].strip() for line in control.splitlines()
                 if line.lstrip().startswith("EMITTER[")]
        missing = [tbl_dir / f"tria_{gas}.nc" for gas in gases
                   if not (tbl_dir / f"tria_{gas}.nc").is_file()]
        if missing:
            parser.error(f"missing lookup tables for {case}: {', '.join(map(str, missing))}")
    tblbase = str(tbl_dir / "tria")
    binary = args.bin.resolve()
    if not binary.is_file():
        parser.error(f"missing executable: {binary}; build with cd src && make -j")
    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    for case in args.cases:
        directory = ROOT / "projects/examples" / case
        for threads in args.threads:
            env = os.environ.copy()
            env.update(OMP_NUM_THREADS=str(threads), LC_ALL="C")
            result = subprocess.run(
                [str(binary), f"{case}.ctl", "obs.tab", "atm.tab", str(output.parent / f"{case}_{threads}.tab"), "TASK", "t", "TBLBASE", tblbase, "TBLFMT", "3"],
                cwd=directory, env=env, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            )
            log = output.parent / f"{case}_{threads}.log"
            log.write_text(result.stdout)
            if result.returncode:
                raise RuntimeError(f"{case}, {threads} threads failed; see {log}")
            match = TIMING.search(result.stdout)
            if not match:
                raise RuntimeError(f"timing line missing from {log}")
            rows.append([case, threads, *match.groups()])
            print(f"{case}: {threads} threads, mean {match.group(1)} s")
    with output.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["case", "threads", "mean_s", "stddev_s", "min_s", "max_s"])
        writer.writerows(rows)
    print(f"Results: {output}")


if __name__ == "__main__":
    main()
