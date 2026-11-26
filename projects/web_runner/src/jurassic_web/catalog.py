"""Discover lookup-table availability from installed or registered datasets."""
from __future__ import annotations
import json
import math
import re
from pathlib import Path

NAME = re.compile(r"^(?P<family>[a-zA-Z][a-zA-Z0-9-]{0,30})_(?P<channel>[0-9]{3,5}\.[0-9]{4})(?:_(?P<gas>[A-Za-z][A-Za-z0-9]{0,15}))?\.(?P<ext>tab|filt)$")


def has_data(path):
    """Basic non-empty-table check, not scientific spectral validation."""
    try:
        with path.open("r", errors="replace") as handle:
            for line in handle:
                line = line.strip()
                if not line or line.startswith("#"): continue
                cols = line.split()
                if len(cols) >= 3:
                    try:
                        if all(math.isfinite(float(x)) for x in cols[:3]): return True
                    except ValueError: pass
                return False
    except (OSError, UnicodeError): pass
    return False


def scan_family(folder, family, source="Built-in"):
    channels = {}
    for path in sorted(folder.glob(f"{family}_*")):
        match = NAME.fullmatch(path.name)
        if not match or match["family"] != family: continue
        ch = match["channel"]
        item = channels.setdefault(ch, {"gases": {}, "filter": False})
        if match["ext"] == "filt": item["filter"] = path.stat().st_size > 0
        if match["ext"] == "tab" and match["gas"]:
            item["gases"][match["gas"]] = {"data": has_data(path), "bytes": path.stat().st_size}
    return {"id": family, "label": family.title(), "source": source, "channels": channels, "prefix": str((folder / family).resolve())}


def datasets(app):
    root = Path(app.config["JURASSIC_ROOT"])
    builtins = root / "tests" / "data"
    output = {}
    if builtins.is_dir():
        for path in builtins.iterdir():
            m = NAME.fullmatch(path.name)
            if m and m["family"] not in output:
                fam = m["family"]
                output[fam] = scan_family(builtins, fam)
    manifests = Path(app.config["INSTANCE_DIR"]) / "datasets"
    if manifests.exists():
        for manifest in sorted(manifests.glob("*/manifest.json")):
            try:
                metadata = json.loads(manifest.read_text())
                family = metadata["family"]
                ds = scan_family(manifest.parent, family, source="Uploaded")
                # ID is independent from the file prefix.
                ds["id"] = manifest.parent.name
                ds["label"] = metadata.get("label", family)
                output[ds["id"]] = ds
            except (OSError, KeyError, ValueError): continue
    return output


def available_gases(ds, channels):
    names = set()
    for ch in channels:
        item = ds["channels"].get(ch, {})
        for gas, desc in item.get("gases", {}).items():
            if desc["data"]: names.add(gas)
    return sorted(names)


def compatible(ds, channels, gases):
    """Return warnings for partial/empty combinations; reject only fully absent pairings."""
    warnings = []
    for ch in channels:
        item = ds["channels"].get(ch)
        if not item: warnings.append(f"Channel {ch} is missing from this dataset."); continue
        if not item["filter"]: warnings.append(f"Missing filter file at {ch} cm⁻¹.")
        present = False
        for gas in gases:
            desc = item["gases"].get(gas)
            if desc and desc["data"]: present = True
            else: warnings.append(f"{gas} at {ch} cm⁻¹: no numerical lookup-table data (or missing file).")
        if not present: warnings.append(f"No usable selected-gas lookup table at {ch} cm⁻¹.")
    return warnings
