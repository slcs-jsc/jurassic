#!/usr/bin/env python3
import argparse
import re


def block_range(lines, pattern):
    idxs = [i for i, l in enumerate(lines) if re.match(pattern, l)]
    return (min(idxs), max(idxs) + 1) if idxs else None


def replace_gases(lines, gases):
    r = block_range(lines, r"^\s*(NG|EMITTER\[\d+\])\s*=")
    block = [f"NG = {len(gases)}\n"]
    block += [f"EMITTER[{i}] = {g}\n" for i, g in enumerate(gases)]
    lines[r[0]:r[1]] = block

    zmin_r = block_range(lines, r"^\s*RETQ_ZMIN\[\d+\]\s*=")
    zmax_r = block_range(lines, r"^\s*RETQ_ZMAX\[\d+\]\s*=")
    if zmin_r and zmax_r:
        zmin = lines[zmin_r[0]].split("=", 1)[1].strip()
        zmax = lines[zmax_r[0]].split("=", 1)[1].strip()
        start, end = min(zmin_r[0], zmax_r[0]), max(zmin_r[1], zmax_r[1])
        block = []
        for i in range(len(gases)):
            block.append(f"RETQ_ZMIN[{i}] = {zmin}\n")
            block.append(f"RETQ_ZMAX[{i}] = {zmax}\n")
        lines[start:end] = block


def replace_channels(lines, nd, nu_min, nu_max):
    r = block_range(lines, r"^\s*(ND|NU\[\d+\])\s*=")
    if nd == 1:
        freqs = [(nu_min + nu_max) / 2]
    else:
        step = (nu_max - nu_min) / (nd - 1)
        freqs = [nu_min + i * step for i in range(nd)]
    block = [f"ND = {nd}\n"]
    block += [f"NU[{i}] = {v:.4f}\n" for i, v in enumerate(freqs)]
    lines[r[0]:r[1]] = block


def main():
    p = argparse.ArgumentParser()
    p.add_argument("template")
    p.add_argument("out")
    p.add_argument("--nd", type=int)
    p.add_argument("--nu-min", type=float, default=900.0)
    p.add_argument("--nu-max", type=float, default=1500.0)
    p.add_argument("--gas-file")
    args = p.parse_args()

    with open(args.template) as f:
        lines = f.readlines()

    if args.gas_file:
        with open(args.gas_file) as f:
            gases = [g.strip() for g in f if g.strip()]
        replace_gases(lines, gases)

    if args.nd is not None:
        replace_channels(lines, args.nd, args.nu_min, args.nu_max)

    with open(args.out, "w") as f:
        f.writelines(lines)


if __name__ == "__main__":
    main()