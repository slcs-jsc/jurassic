import re
from pathlib import Path
import statistics as st

# Matches the file naming produced by base.sh's jr_run():
#   <label>.t<threads>.<group>.b<batch>.rep<rep>.csv
# e.g. "noise.t1.MEM_DP.b240.rep3.csv", "base.t24.CACHE.b240.rep1.csv"
RUN_FILE_RE = re.compile(
    r"^(?P<label>[A-Za-z0-9_]+)\."
    r"t(?P<threads>\d+)\."
    r"(?P<group>[A-Za-z0-9_]+)\."
    r"b(?P<batch>\d+)\."
    r"rep(?P<rep>\d+)"
    r"\.(?P<ext>csv|txt)$"
)

PERMISSION_RE = re.compile(r"Setup of event (\S+) on CPU \d+ failed: Permission denied")
TIMER_RE = re.compile(r"^(TIMER_\w+)\s*=\s*([\d.eE+-]+)\s*s", re.MULTILINE)
RUNTIME_RE = re.compile(
    r"RUNTIME:\s*execution=\s*(\S+)\s*\|\s*batch_size=\s*(\d+)\s*\|\s*"
    r"mean=\s*([\d.eE+-]+)\s*s\s*\|\s*stddev=\s*([\d.eE+-]+)\s*s\s*\|\s*"
    r"min=\s*([\d.eE+-]+)\s*s\s*\|\s*max=\s*([\d.eE+-]+)\s*s"
)

def _to_number(s: str):
    try:
        return float(s)
    except ValueError:
        return s

def parse_likwid_profile(path: Path) -> dict:
    """
    Parse one likwid-perfctr CSV (marker mode, -m).
 
    Returns:
        {
          "cpu_info": {...},
          "regions": {
              "<region name, e.g. 'formod'>": {
                  "region_info": {metric_name: [per-thread values]},
                  "tables": {"<group>:<kind>": {metric_name: [values]}},
              },
              ...
          },
        }
    """
    lines = Path(path).read_text().splitlines()
    cpu_info = {}
    regions: dict[str, dict] = {}
    current_region = ""
    i = 0
    while i < len(lines):
        row = [c.strip() for c in lines[i].split(",")]
        if not row or not row[0]:
            i += 1
            continue

        if row[0] == "STRUCT":
            _, name, count = row[0], row[1], int(row[2])
            for j in range(count):
                key, val = lines[i + 1 + j].split(",")[:2]
                cpu_info[key.strip().rstrip(":")] = val.strip()
            i += 1 + count

        elif row[0] == "TABLE":
            fields = [c for c in row]
            while fields and fields[-1] == "":
                fields.pop()

            region_field = next((f for f in fields if f.startswith("Region ")), None)
            if region_field is not None:
                current_region = region_field[len("Region "):].strip()
                fields.remove(region_field)
            regions.setdefault(current_region, {"region_info": {}, "tables": {}})

            count = int(fields[-1])
            group = fields[-2]
            kind = fields[-3]
            i += 1

            # Region Info block: the "Region Info" line IS the header row
            # (label + per-thread columns, same orientation as the main
            # metric table below), followed by exactly two fixed data rows:
            # "RDTSC Runtime [s]" and "call count".
            peek = [c.strip() for c in lines[i].split(",")]
            while peek and peek[-1] == "":
                peek.pop()
            if peek and peek[0] == "Region Info":
                thread_cols = peek[1:]
                for j in range(2):
                    cells = [c.strip() for c in lines[i + 1 + j].split(",")]
                    name = cells[0]
                    values = cells[1:1 + len(thread_cols)]
                    regions[current_region]["region_info"][name] = [
                        _to_number(v) for v in values
                    ]
                i += 1 + 2

            header = [c.strip() for c in lines[i].split(",") if c.strip()]
            thread_cols = header[1:]
            data = {}
            for j in range(count):
                cells = [c.strip() for c in lines[i + 1 + j].split(",")]
                metric_name, values = cells[0], cells[1:1 + len(thread_cols)]
                data[metric_name] = [_to_number(v) for v in values]
            regions[current_region]["tables"][f"{group}:{kind}"] = data
            i += 1 + count
        else:
            i += 1

    return {"cpu_info": cpu_info, "regions": regions}

def parse_formod_log(path: Path) -> dict:
    text = Path(path).read_text()
    permission_errors = sorted(set(PERMISSION_RE.findall(text)))
 
    batch = None
    m = RUNTIME_RE.search(text)
    if m:
        batch = {
            "execution": m.group(1), "batch_size": int(m.group(2)),
            "mean_s": float(m.group(3)), "stddev_s": float(m.group(4)),
            "min_s": float(m.group(5)), "max_s": float(m.group(6)),
        }
 
    timers = {name: float(val) for name, val in TIMER_RE.findall(text)}
 
    return {"permission_errors": permission_errors, "batch": batch, "timers": timers}
 
def parse_run_dir(run_dir: Path) -> list:
    """
    Parse every <label>.t<N>.<group>.b<N>.rep<N>.csv in run_dir/out/
    (as produced by common.sh's jr_run/jr_finish).
 
    Each returned entry corresponds to one CSV file and carries:
        label, threads, group, batch, rep, group_recognized,
        cpu_info, regions (see parse_likwid_profile), and, if a matching
        .txt exists, permission_errors / batch / timers from parse_formod_log.
    """
    run_dir = Path(run_dir)
    out_dir = run_dir / "out"
    search_dir = out_dir if out_dir.is_dir() else run_dir
 
    available_groups_file = run_dir / "likwid_available_groups.txt"
    available_groups = available_groups_file.read_text() if available_groups_file.exists() else ""
 
    configs = []
    for csv_path in sorted(search_dir.glob("*.csv")):
        match = RUN_FILE_RE.match(csv_path.name)
        if not match:
            continue
 
        label = match.group("label")
        threads = int(match.group("threads"))
        group = match.group("group")
        batch_size = int(match.group("batch"))
        rep = int(match.group("rep"))
        txt_path = csv_path.with_suffix(".txt")
 
        csv_data = parse_likwid_profile(csv_path)
 
        entry = {
            "label": label,
            "threads": threads,
            "group": group,
            "batch_size": batch_size,
            "rep": rep,
            "group_recognized": bool(
                re.search(rf"^\s*{re.escape(group)}\b", available_groups, re.MULTILINE)
            ),
            "cpu_info": csv_data["cpu_info"],
            "regions": csv_data["regions"],
        }
        if txt_path.exists():
            entry.update(parse_formod_log(txt_path))
        configs.append(entry)
 
    return configs

def per_call(raw, call_count):
    if raw is None or call_count in (None, 0) or isinstance(raw, list):
        return None
    return raw / call_count

def get_stats(values):
    mean = st.mean(values)
    median = st.median(values)
    stdev = st.pstdev(values)
    if mean == 0 : 
        cv = float("nan")
    else: 
        cv = stdev / mean
    return mean, median, stdev, cv

CAS_LINE_SIZE_BYTES = 64
def get_cas_total_gbytes(entry, region, prefix):
    """CAS_COUNT_RD/WR appear once per MBOX channel, not as a single named row.
    Need to sum over channels and convert to GBytes (1 CAS = one 64-byte transfer)."""
    tables = entry.get("regions", {}).get(region, {}).get("tables", {})
    group_prefix = f"{entry['group']}:"
    total_cas = 0
    found = False
    for key, table in tables.items():
        if not key.startswith(group_prefix):
            continue
        for name, values in table.items():
            if not name.startswith(prefix):
                continue
            if isinstance(values, list):
                numeric = [v for v in values if isinstance(v, (int, float))]
                if not numeric:
                    continue
                total_cas += numeric[0]
                found = True
            elif isinstance(values, (int, float)):
                total_cas += values
                found = True
    return total_cas * CAS_LINE_SIZE_BYTES / 1e9 if found else None

def get_value(entry, region, metric):
    if metric in ("CAS_COUNT_RD", "CAS_COUNT_WR"):
        return get_cas_total_gbytes(entry, region, metric)
    return get_metric(entry, region, metric)
 
def get_metric(entry: dict, region: str, metric: str, stat: bool = False):
    """
    Searches any table belonging to `entry["group"]` in the given region
    """
    tables = entry.get("regions", {}).get(region, {}).get("tables", {})
    prefix = f"{entry['group']}:"

    for key, table in tables.items():
        if not key.startswith(prefix):
            continue
        kind = key[len(prefix):].strip()
        is_stat_table = kind.upper().endswith("STAT")
        if is_stat_table != stat:
            continue
        if metric in table:
            values = table[metric]
            if stat:
                keys = ["sum", "min", "max", "avg"]
                return dict(zip(keys, values))
            numeric = [v for v in values if isinstance(v, (int, float))]
            return numeric[0] if len(numeric) == 1 else sum(numeric)
    return None

def normalize(raw, metric, entry, call_count):
    m = metric.lower()
    if "bandwidth" in m or "mflop/s" in m: 
        # already a rate, socket-wide
        return raw
    if "volume" in m or "energy" in m:
        # socket-wide total ÷ scenes
        return raw / entry["batch_size"]
    return per_call(raw, call_count)

def get_summed_metric(entry, region, prefix):
    """
    Helper function to sum up matching rows
    """
    tables = entry.get("regions", {}).get(region, {}).get("tables", {})
    key = next((k for k in tables if k.startswith(f"{entry['group']}:")), None)
    if not key:
        return None
    total = 0
    found = False
    for name, values in tables[key].items():
        if name.startswith(prefix):
            total += values[0] if not isinstance(values, list) else sum(v for v in values if isinstance(v, (int, float)))
            found = True
    return total if found else None
 
def get_call_count(entry: dict, region: str):
    """Total call count recorded for `region` (summed across threads)."""
    info = entry.get("regions", {}).get(region, {}).get("region_info", {})
    for key in info:
        if key.lower().startswith("call count"):
            vals = info[key]
            return sum(v for v in vals if isinstance(v, (int, float)))
    return None
 
def get_region_runtime(entry: dict, region: str):
    """Total RDTSC runtime recorded for `region` (summed across threads)."""
    info = entry.get("regions", {}).get(region, {}).get("region_info", {})
    for key in info:
        if "runtime" in key.lower():
            vals = info[key]
            return sum(v for v in vals if isinstance(v, (int, float)))
    return None

def collect(entries, region, metric, warmup):
    """Median + CV of one metric across the kept repetitions."""
    entries = sorted(entries, key=lambda e: e["rep"])[warmup:]
    vals = []
    for e in entries:
        raw = get_value(e, region, metric)
        if raw is None:
            continue
        v = normalize(raw, metric, e, get_call_count(e, region))
        if v is not None:
            vals.append(v)
    if not vals:
        return None, None, None, float("nan"), 0, []
    
    return *get_stats(vals), len(vals), vals
 
def collect_runtime(entries, warmup):
    entries = sorted(entries, key=lambda e: e["rep"])[warmup:]
    vals = []
    for e in entries:
        b = e.get("batch")
        if b:
            vals.append(b["mean_s"])
    if not vals:
        return None, None, None, float("nan"), 0, []
    return *get_stats(vals), len(vals), vals

def collect_intensity(entries, region, warmup):
    """
    Median + CV of arithmetic intensity across reps
    computed as (summed DP flops) / (summed memory bandwidth)
    """
    entries = sorted(entries, key=lambda e: e["rep"])[warmup:]
    vals = []
    for e in entries:
        flops = get_value(e, region, "DP [MFLOP/s]")
        bw = get_value(e, region, "Memory bandwidth [MBytes/s]")
        if flops is not None and bw:
            vals.append(flops / bw)
    if not vals:
        return None, None, None, float("nan"), 0, []
    return *get_stats(vals), len(vals), vals
 
def summarize(configs: list) -> str:
    if not configs:
        return "(no benchmark configs found)"
 
    perm_errors = sorted({e for c in configs for e in c.get("permission_errors", [])})
    lines = []
    if perm_errors:
        lines.append(
            f"NOTE: counters unavailable on this node (permission denied): "
            f"{', '.join(perm_errors)}. Related metrics read 0 and are not "
            f"real measurements."
        )
    unrecognized = sorted({c["group"] for c in configs if not c["group_recognized"]})
    if unrecognized:
        lines.append(f"NOTE: requested group(s) not recognized on this node: "
                      f"{', '.join(unrecognized)}.")
 
    lines.append(
        f"{'label':>10} | {'threads':>7} | {'group':>8} | {'batch':>6} | "
        f"{'rep':>4} | {'analysis_s':>10} | regions"
    )
    for c in sorted(configs, key=lambda c: (c["label"], c["threads"], c["group"], c["rep"])):
        timers = c.get("timers", {})
        analysis = timers.get("TIMER_GROUP_ANALYSIS", "-")
        region_names = ", ".join(sorted(c.get("regions", {}).keys())) or "-"
        lines.append(
            f"{c['label']:>10} | {c['threads']:>7} | {c['group']:>8} | "
            f"{c['batch_size']:>6} | {c['rep']:>4} | {analysis!s:>10} | {region_names}"
        )
 
    return "\n".join(lines)