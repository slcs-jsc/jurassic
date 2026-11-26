"""Constrained control-file parsing and generation. No arbitrary path directives."""
from __future__ import annotations
import math
import re
from pathlib import Path
from .settings import SETTINGS

KEY = re.compile(r"^([A-Za-z][A-Za-z0-9_]*(?:\[(?:\d+|\*)\])?)\s*=\s*(.*?)\s*$")
MANAGED = re.compile(r"^(?:NG|ND|TBLBASE|EMITTER\[\d+\]|NU\[\d+\]|RETQ_Z(?:MIN|MAX)\[\d+\]|WRITE_BBT|RAYDS|RAYDZ|REFRAC)$")
ALLOWED_IMPORT = re.compile(r"^(?:NG|ND|TBLBASE|EMITTER\[\d+\]|NU\[\d+\]|WRITE_BBT|RAYDS|RAYDZ|REFRAC)$")
GAS = re.compile(r"^[A-Za-z][A-Za-z0-9]{0,15}$")


def parse_ctl(text, imported=False):
    if len(text) > 100_000 or "\x00" in text: raise ValueError("Invalid or oversized control file")
    data = {}
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"): continue
        m = KEY.fullmatch(line)
        if not m:
            if imported: raise ValueError("Invalid control-file syntax")
            continue
        k,v = m.groups()
        if imported and not ALLOWED_IMPORT.fullmatch(k):
            continue  # Never import filesystem paths or unreviewed options
        data.setdefault(k,v)
    return data


def import_form(text):
    x = parse_ctl(text, imported=True)
    if not x: raise ValueError("No supported configuration parameters found")
    def numbered(prefix, n):
        count = min(SETTINGS["max_channels"] if prefix == "NU" else SETTINGS["max_gases"], max(0, int(x.get(n, "0"))))
        return [x[f"{prefix}[{i}]"] for i in range(count) if f"{prefix}[{i}]" in x]
    result = {"gases": numbered("EMITTER", "NG"), "channels": numbered("NU", "ND")}
    result["output"] = "bbt" if x.get("WRITE_BBT", "0") == "1" else "radiance"
    for key in ("RAYDS", "RAYDZ", "REFRAC"):
        if key in x: result[key.lower()] = x[key]
    return result


def normalized_channels(channels):
    result = []
    for value in channels:
        try: f = float(value)
        except (ValueError, TypeError): raise ValueError("Channel wavenumbers must be numeric")
        if not math.isfinite(f) or not 100 <= f <= 5000: raise ValueError("Channel wavenumber out of allowed range")
        v = f"{f:.4f}"
        if v not in result: result.append(v)
    return result


def as_float(v, low, high, label):
    try: f = float(v)
    except (ValueError, TypeError): raise ValueError(f"{label} must be a number")
    if not math.isfinite(f) or f < low or f > high:
        raise ValueError(f"{label} must be between {low} and {high}")
    return f


def normalize(raw, ds):
    if not isinstance(raw, dict): raise ValueError("Invalid simulation request")
    geo = str(raw.get("geometry", "")).lower()
    if geo not in SETTINGS["geometries"]: raise ValueError("Choose Nadir, Limb or Zenith")
    channels = normalized_channels(raw.get("channels", []))
    gases = list(dict.fromkeys(raw.get("gases", [])))
    if not channels or len(channels) > SETTINGS["max_channels"]: raise ValueError("Select 1–16 spectral channels")
    if not gases or len(gases) > SETTINGS["max_gases"] or any(not GAS.fullmatch(g) for g in gases):
        raise ValueError("Select 1–16 valid gases")
    if any(ch not in ds["channels"] for ch in channels): raise ValueError("A selected channel is not in this dataset")
    if any(not any(gas in ds["channels"][ch]["gases"] for ch in channels) for gas in gases):
        raise ValueError("A selected gas is not listed for the selected channels")
    if any(not ds["channels"][ch]["filter"] for ch in channels): raise ValueError("Missing channel filter file")
    if any(not any(ds["channels"][ch]["gases"].get(g, {}).get("data") for g in gases) for ch in channels):
        raise ValueError("At least one usable gas lookup table is needed for each channel")
    output = raw.get("output", "radiance")
    if output not in {"bbt", "radiance"}: raise ValueError("Invalid output quantity")
    rds = as_float(raw.get("rayds", 10), 0.1, 100, "RAYDS")
    rdz = as_float(raw.get("raydz", 0.1), 0.01, 5, "RAYDZ")
    refrac = raw.get("refrac", True)
    if not isinstance(refrac, bool): raise ValueError("Invalid refraction setting")
    kernel = raw.get("kernel", False)
    if not isinstance(kernel, bool): raise ValueError("Invalid kernel setting")
    if kernel and (len(channels)*len(gases) > 20): raise ValueError("Kernel configuration exceeds GUI resource limits")
    return {"geometry": geo, "dataset": ds["id"], "gases": gases, "channels": channels,
            "output": output, "rayds": rds, "raydz": rdz, "refrac": refrac, "kernel": kernel}


def make_control(root, parameters, ds):
    geometry = parameters["geometry"]
    example = Path(root) / "projects" / "examples" / geometry / f"{geometry}.ctl"
    if not example.is_file(): raise ValueError(f"Missing upstream example control file: {geometry}")
    raw = example.read_text()
    kept = []
    for line in raw.splitlines():
        match = KEY.fullmatch(line.strip())
        if match and MANAGED.fullmatch(match[1]): continue
        kept.append(line)
    dynamic = ["", "# Web Runner generated values (validated)", f"TBLBASE = {ds['prefix']}",
               f"NG = {len(parameters['gases'])}"]
    dynamic += [f"EMITTER[{i}] = {gas}" for i,gas in enumerate(parameters["gases"])]
    dynamic += [f"RETQ_ZMIN[{i}] = -100\nRETQ_ZMAX[{i}] = 88" for i,_ in enumerate(parameters["gases"])]
    dynamic += [f"ND = {len(parameters['channels'])}"]
    dynamic += [f"NU[{i}] = {channel}" for i,channel in enumerate(parameters["channels"])]
    dynamic += [f"WRITE_BBT = {1 if parameters['output']=='bbt' else 0}",
                f"RAYDS = {parameters['rayds']:g}", f"RAYDZ = {parameters['raydz']:g}",
                f"REFRAC = {int(parameters['refrac'])}"]
    return "\n".join(kept + dynamic) + "\n"
