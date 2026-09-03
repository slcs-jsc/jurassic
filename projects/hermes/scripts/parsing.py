import re
from pathlib import Path

CSV_RE = re.compile(r"log\.omp(\d+)\.([A-Za-z0-9_]+)\.csv$")
PERMISSION_RE = re.compile(r"Setup of event (\S+) on CPU \d+ failed: Permission denied")
TIMER_RE = re.compile(r"^(TIMER_\w+)\s*=\s*([\d.eE+-]+)\s*s", re.MULTILINE)
RUNTIME_RE =  re.compile(
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
    lines = Path(path).read_text().splitlines()
    cpu_info, tables = {}, {}
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
            # Trim trailing empty padding fields
            # likwid-marker-Mode inserts an extra "Region <name>" field after TABLE so the field count varies
            fields = [c for c in row]
            while fields and fields[-1] == "":
                fields.pop()
            count = int(fields[-1])
            group = fields[-2]
            kind = fields[-3]
            i += 1

            # likwid-marker-Mode insert a one-off "Region Info" mini-block
            peek = [c.strip() for c in lines[i].split(",")]
            if peek and peek[0] == "Region Info":
                i += 3 

            header = [c.strip() for c in lines[i].split(",") if c.strip()]
            thread_cols = header[1:]
            data = {}
            for j in range(count):
                cells = [c.strip() for c in lines[i + 1 + j].split(",")]
                metric_name, values = cells[0], cells[1:1 + len(thread_cols)]
                data[metric_name] = [_to_number(v) for v in values]
            tables[f"{group}:{kind}"] = data
            i += 1 + count
        else:
            i += 1
    return {"cpu_info": cpu_info, "tables": tables}

def parse_formod_log(path:Path) -> dict:

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
    run_dir = Path(run_dir)
    available_groups_file = run_dir / "likwid_available_groups.txt"
    available_groups = available_groups_file.read_text() if available_groups_file.exists() else ""
 
    configs = []
    for csv_path in sorted(run_dir.glob("log.omp*.csv")):
        match = CSV_RE.search(csv_path.name)
        if not match:
            continue
        omp_threads, group = int(match.group(1)), match.group(2)
        txt_path = csv_path.with_suffix(".txt")
 
        csv_data = parse_likwid_profile(csv_path)
        metric_table = csv_data["tables"].get(f"{group}:Group 1 Metric", {})
        metrics = {k: v[0] if len(v) == 1 else v for k, v in metric_table.items()}

        stat_table = csv_data["tables"].get(f"{group}:Group 1 Metric STAT", {})
        metrics_stat = {
            metric_name: dict(zip(["sum", "min", "max", "avg"], values))
            for metric_name, values in stat_table.items()
        }

        entry = {
            "omp_threads": omp_threads,
            "group": group,
            "group_recognized": bool(re.search(rf"^\s*{re.escape(group)}\b", available_groups, re.MULTILINE)),
            "cpu_info": csv_data["cpu_info"],
            "metrics": metrics,
            "metrics_stat": metrics_stat,
        }
        if txt_path.exists():
            entry.update(parse_formod_log(txt_path))
        configs.append(entry)

    return configs
 
def summarize(configs: list) -> str: 

    if not configs:
        return "(no benchmark configs found)"
 
    perm_errors = sorted({e for c in configs for e in c.get("permission_errors", [])})
    lines = []
    if perm_errors:
        lines.append(
            f"NOTE: counters unavailable on this node (permission denied): {', '.join(perm_errors)}. "
            "Related metrics (e.g. Clock [MHz], Runtime unhalted) read 0 and are not real measurements."
        )
    unrecognized = sorted({c["group"] for c in configs if not c["group_recognized"]})
    if unrecognized:
        lines.append(f"NOTE: requested group(s) not recognized on this node: {', '.join(unrecognized)}.")
 
    lines.append(f"{'threads':>7} | {'group':>8} | {'analysis_s':>10} | {'input_s':>8} | {'total_s':>8} | key metrics")
    for c in sorted(configs, key=lambda c: (c["omp_threads"], c["group"])):
        timers = c.get("timers", {})
        analysis = timers.get("TIMER_GROUP_ANALYSIS", "-")
        inp = timers.get("TIMER_GROUP_INPUT", "-")
        total = timers.get("TIMER_TOTAL", "-")
        key_metrics = ", ".join(f"{k}={v}" for k, v in c["metrics"].items() if k != "Runtime (RDTSC) [s]")
        lines.append(f"{c['omp_threads']:>7} | {c['group']:>8} | {analysis!s:>10} | {inp!s:>8} | {total!s:>8} | {key_metrics}")
 
    return "\n".join(lines)