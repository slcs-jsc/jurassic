#!/usr/bin/env python3
"""Write a control file with ND channels taken from a channel list and the
emitters of a gas file. Prints the number of active (channel, gas) pairs."""
import argparse
import re
import sys


def block_range(lines, pattern):
    idxs = [i for i, l in enumerate(lines) if re.match(pattern, l)]
    return (min(idxs), max(idxs) + 1) if idxs else None


def read_gases(path):
    with open(path) as f:
        return [l.strip() for l in f if l.strip() and not l.lstrip().startswith("#")]


def read_channels(path):
    """Rows of the channel list: (wavenumber, set of active gases)."""
    rows = []
    with open(path) as f:
        next(f)
        for line in f:
            if line.strip():
                nu, _, gases = line.rstrip("\n").split("\t")
                rows.append((int(nu), set(gases.split())))
    return rows


def pick_channels(rows, nd):
    """nd evenly spaced rows, first and last included."""
    if not 1 <= nd <= len(rows):
        sys.exit(f"ND={nd} not possible with {len(rows)} channels in the list")
    if nd == 1:
        return [rows[len(rows) // 2]]
    return [rows[round(i * (len(rows) - 1) / (nd - 1))] for i in range(nd)]


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


def replace_channels(lines, nus):
    r = block_range(lines, r"^\s*(ND|NU\[\d+\])\s*=")
    block = [f"ND = {len(nus)}\n"]
    block += [f"NU[{i}] = {nu}\n" for i, nu in enumerate(nus)]
    lines[r[0]:r[1]] = block


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("template")
    p.add_argument("out")
    p.add_argument("--channels", required=True, help="channel list (NU, NGAS, GASES)")
    p.add_argument("--nd", type=int, required=True)
    p.add_argument("--gas-file", required=True)
    args = p.parse_args()

    with open(args.template) as f:
        lines = f.readlines()

    gases = read_gases(args.gas_file)
    channels = pick_channels(read_channels(args.channels), args.nd)
    replace_gases(lines, gases)
    replace_channels(lines, [nu for nu, _ in channels])

    with open(args.out, "w") as f:
        f.writelines(lines)
    print(sum(len(active & set(gases)) for _, active in channels))


if __name__ == "__main__":
    main()
