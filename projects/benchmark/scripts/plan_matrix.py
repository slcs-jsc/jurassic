#!/usr/bin/env python3
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import csv
import signal
import sys

ROOT = Path(__file__).resolve().parents[1]
CONFIGS = ROOT / "configs"

GEOMETRIES = ["zenith", "nadir", "limb"]
CPU_THREADS = [1, 2, 4, 8, 12]
GPU_BATCHES = [1, 8, 64, 256]
CHANNEL_COUNTS = [int(x.strip()) for x in (CONFIGS / "channel_counts.txt").read_text().splitlines() if x.strip()]


@dataclass(frozen=True)
class Family:
    name: str
    geometries: tuple[str, ...]
    band_names: tuple[str, ...]
    channel_counts: tuple[int, ...]
    gas_sets: tuple[str, ...]
    cpu_threads: tuple[int, ...]
    gpu_batches: tuple[int, ...]
    ray_mode: str
    gpu_only: bool = False
    cpu_only: bool = False


RAY_SETS: dict[str, dict[str, int]] = {}
with (CONFIGS / "rays.tsv").open() as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        RAY_SETS.setdefault(row["geometry"], {})[row["ray_set"]] = int(row["nr"])


FAMILIES = {
    "geometry_baseline": Family(
        name="geometry_baseline",
        geometries=tuple(GEOMETRIES),
        band_names=("band_0900_1500",),
        channel_counts=(32,),
        gas_sets=("core",),
        cpu_threads=(1, 12),
        gpu_batches=(1, 8, 64),
        ray_mode="medium",
    ),
    "channel_scaling": Family(
        name="channel_scaling",
        geometries=tuple(GEOMETRIES),
        band_names=("band_0900_1500",),
        channel_counts=tuple(CHANNEL_COUNTS),
        gas_sets=("core",),
        cpu_threads=(1, 12),
        gpu_batches=(1, 8, 64, 256),
        ray_mode="medium",
    ),
    "gas_scaling": Family(
        name="gas_scaling",
        geometries=tuple(GEOMETRIES),
        band_names=("band_0900_1500",),
        channel_counts=(64, 256),
        gas_sets=("core", "priority_mid", "priority_full"),
        cpu_threads=(1, 12),
        gpu_batches=(1, 8, 64),
        ray_mode="medium",
    ),
    "spectral_band_scaling": Family(
        name="spectral_band_scaling",
        geometries=tuple(GEOMETRIES),
        band_names=("band_0500_0900", "band_0900_1500", "band_1500_2200", "band_2200_3000"),
        channel_counts=(64, 256),
        gas_sets=("core", "priority_full"),
        cpu_threads=(1,),
        gpu_batches=(1, 8, 64),
        ray_mode="medium",
    ),
    "cpu_strong_scaling": Family(
        name="cpu_strong_scaling",
        geometries=tuple(GEOMETRIES),
        band_names=("band_0900_1500",),
        channel_counts=(64,),
        gas_sets=("core",),
        cpu_threads=tuple(CPU_THREADS),
        gpu_batches=(),
        ray_mode="medium",
        cpu_only=True,
    ),
    "gpu_throughput_scaling": Family(
        name="gpu_throughput_scaling",
        geometries=("zenith", "nadir", "limb"),
        band_names=("band_0900_1500",),
        channel_counts=(1, 16, 64, 256),
        gas_sets=("core", "priority_mid"),
        cpu_threads=(),
        gpu_batches=tuple(GPU_BATCHES),
        ray_mode="custom_gpu",
        gpu_only=True,
    ),
    "ray_scaling": Family(
        name="ray_scaling",
        geometries=tuple(GEOMETRIES),
        band_names=("band_0900_1500",),
        channel_counts=(64,),
        gas_sets=("core",),
        cpu_threads=(1, 12),
        gpu_batches=(1, 8, 64),
        ray_mode="all",
    ),
}


def resolve_ray_sets(geometry: str, mode: str):
    if mode == "medium":
        return [("medium", RAY_SETS[geometry]["medium"])]
    if mode == "all":
        return [(name, RAY_SETS[geometry][name]) for name in ("small", "medium", "large")]
    if mode == "custom_gpu":
        custom = {
            "zenith": "medium",
            "nadir": "small",
            "limb": "medium",
        }[geometry]
        return [(custom, RAY_SETS[geometry][custom])]
    raise ValueError(mode)


def emit_family(family: Family):
    rows = []
    for geometry in family.geometries:
        for ray_set, nr in resolve_ray_sets(geometry, family.ray_mode):
            for band in family.band_names:
                for nd in family.channel_counts:
                    for gas_set in family.gas_sets:
                        if not family.gpu_only:
                            for omp in family.cpu_threads:
                                rows.append((family.name, geometry, ray_set, nr, band, nd, gas_set, "cpu", omp, ""))
                        if not family.cpu_only:
                            for batch in family.gpu_batches:
                                rows.append((family.name, geometry, ray_set, nr, band, nd, gas_set, "gpu", "", batch))
    return rows


def main(argv: list[str]) -> int:
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    selected = argv[1:] if len(argv) > 1 else list(FAMILIES)
    for name in selected:
        if name not in FAMILIES:
            print(f"Unknown family: {name}", file=sys.stderr)
            return 2

    print("family\tgeometry\tray_set\tnr\tband\tnd\tgas_set\ttarget\tomp_threads\tbatch_size")
    for name in selected:
        for row in emit_family(FAMILIES[name]):
            print("\t".join(str(x) for x in row))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
