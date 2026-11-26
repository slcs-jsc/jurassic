"""Read JURASSIC ASCII result files and prepare chart/table data."""
from __future__ import annotations
import re
from pathlib import Path
import numpy as np

COLUMN = re.compile(r"^#\s*\$(\d+)\s*=\s*(.+)$")


def read_table(path, limit=10000):
    path = Path(path)
    if not path.exists(): raise ValueError(f"Missing result file: {path.name}")
    columns = []
    with path.open("r", errors="replace") as f:
        for line in f:
            m = COLUMN.match(line)
            if m: columns.append(m[2].strip())
    try:
        from itertools import islice
        with path.open("r", errors="replace") as f:
            lines = (line for line in f if line.strip() and not line.lstrip().startswith("#"))
            data = np.loadtxt(islice(lines, limit), ndmin=2)
    except (OSError, ValueError) as exc:
        raise ValueError(f"Could not read {path.name}: {exc}") from exc
    if data.size == 0 or data.ndim != 2: raise ValueError(f"No numeric rows in {path.name}")
    if not np.isfinite(data).all(): raise ValueError(f"Non-finite values in {path.name}")
    return data, columns


def summary(path, channels):
    data, headers = read_table(path)
    expected = 10 + 2 * len(channels)
    if data.shape[1] != expected:
        raise ValueError(f"Unexpected number of columns: {data.shape[1]} (expected {expected})")
    return {"rows": int(data.shape[0]), "columns": int(data.shape[1]), "headers": headers,
            "value_min": float(data[:,10:10+len(channels)].min()),
            "value_max": float(data[:,10:10+len(channels)].max()),
            "transmittance_min": float(data[:,10+len(channels):].min()),
            "transmittance_max": float(data[:,10+len(channels):].max())}


def chart(path, geometry, channels, metric="radiance", max_rows=1000):
    data, headers = read_table(path, limit=max_rows)
    nd = len(channels)
    if data.shape[1] != 10 + 2*nd: raise ValueError("Column count does not match selected channels")
    if metric not in ("radiance", "transmittance"): raise ValueError("Unknown chart quantity")
    xlabel, xidx = {"nadir": ("View-point latitude (°)", 6),
                    "limb": ("Tangent altitude (km)", 7),
                    "zenith": ("View-point latitude (°)", 6)}[geometry]
    ylabel = "Transmittance (–)" if metric == "transmittance" else (headers[10] if len(headers)>10 else "Radiance / brightness temperature")
    idx = 10+nd if metric == "transmittance" else 10
    return {"x": data[:,xidx].tolist(), "xlabel": xlabel, "ylabel": ylabel,
            "traces": [{"name": f"{ch} cm⁻¹", "y": data[:,idx+i].tolist()} for i,ch in enumerate(channels)],
            "headers": headers, "metric": metric}


def rows_preview(path, max_rows=40):
    data, headers = read_table(path, limit=max_rows)
    return {"headers": headers or [f"Column {i+1}" for i in range(data.shape[1])],
            "rows": [[float(v) for v in row] for row in data], "preview_rows": len(data)}


def atmospheric_chart(path, gases):
    data, headers = read_table(path, limit=2000)
    if data.shape[1] < 6: raise ValueError("Invalid atmospheric profile")
    traces = [
        {"name":"Temperature (K)","x":data[:,5].tolist(),"y":data[:,1].tolist(),"type":"temperature"},
        {"name":"Pressure (hPa)","x":data[:,4].tolist(),"y":data[:,1].tolist(),"type":"pressure"},
    ]
    for i, gas in enumerate(gases):
        if 6+i < data.shape[1]:
            traces.append({"name":f"{gas} mixing ratio", "x":data[:,6+i].tolist(), "y":data[:,1].tolist(),"type":"gas"})
    return {"traces": traces, "xlabel":"Value", "ylabel":"Altitude (km)", "headers":headers}


def kernel_heatmap(path, quantity="TEMPERATURE", max_rows=110, max_altitudes=150):
    """Create a sampled heatmap from the ASCII kernel matrix's 13-column format."""
    from collections import defaultdict
    values = defaultdict(dict)
    altitudes = set()
    with Path(path).open("r", errors="replace") as stream:
        for line in stream:
            if not line.strip() or line.startswith("#"): continue
            fields = line.split()
            if len(fields) != 13 or fields[7] != quantity: continue
            try:
                measurement = int(fields[0]); altitude = float(fields[9]); sensitivity = float(fields[12])
            except ValueError: continue
            if measurement >= max_rows or not np.isfinite(sensitivity): continue
            altitudes.add(altitude)
            values[measurement][altitude] = sensitivity
    heights = sorted(altitudes)[:max_altitudes]
    indexes = sorted(values)
    if not heights or not indexes: raise ValueError("No kernel entries for " + quantity)
    return {"x": heights, "y": indexes,
            "z": [[values[i].get(z) for z in heights] for i in indexes],
            "xlabel":"Altitude (km)","ylabel":"Measurement index", "quantity":quantity}
