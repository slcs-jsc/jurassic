"""Local bounded background execution. Replace with external queue on public VM."""
from __future__ import annotations
import json
import os
import signal
import subprocess
import time
import traceback
import uuid
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from .catalog import datasets
from .configuration import make_control
from .parsers import summary
from .plotting import render_chart, render_atmosphere, render_kernel
from .storage import now


def save_atmosphere(data: bytes, path: Path, gases):
    if not data or len(data) > 3 * 1024 * 1024 or b"\x00" in data:
        raise ValueError("Atmospheric file must be a text .tab file below 3 MB")
    try: text = data.decode("utf-8")
    except UnicodeDecodeError: raise ValueError("Atmospheric file must be UTF-8 text")
    headers = [x[1:].strip() for x in text.splitlines() if x.startswith("#")]
    if len(headers) < 6 or not all(any(g.lower() in h.lower() for h in headers) for g in gases):
        raise ValueError("Atmospheric header must name all selected gases")
    numeric = [x.strip() for x in text.splitlines() if x.strip() and not x.startswith("#")]
    if not numeric or len(numeric) > 2000: raise ValueError("Atmosphere must contain 1–2000 rows")
    required = 7 + len(gases)  # time, z, lon, lat, p, t, gases, >=1 window
    for line in numeric:
        fields = line.split()
        if len(fields) < required: raise ValueError("Atmospheric columns do not match selected gases")
        try:
            import math
            if not all(math.isfinite(float(v)) for v in fields): raise ValueError("Atmosphere contains nonfinite values")
        except ValueError as exc: raise ValueError("Atmospheric row contains invalid numbers") from exc
    path.write_bytes(data)


class RunQueue:
    def __init__(self, app):
        self.app = app
        self.pool = ThreadPoolExecutor(max_workers=1, thread_name_prefix="jurassic-run")

    def submit(self, params, atmosphere=None):
        repo = self.app.extensions["jurassic_repo"]
        if repo.active_count() >= self.app.config["MAX_ACTIVE_RUNS"]:
            raise ValueError("Run queue is full; please try again after a run finishes")
        uid = uuid.uuid4().hex[:16]
        run_dir = Path(self.app.config["INSTANCE_DIR"]) / "runs" / uid
        run_dir.mkdir(parents=True, exist_ok=False)
        if atmosphere is not None:
            save_atmosphere(atmosphere, run_dir / "atm.tab", params["gases"])
            params = {**params, "uploaded_atmosphere": True}
        else:
            params = {**params, "uploaded_atmosphere": False}
        repo.create(uid, params)
        try: self.pool.submit(self.execute, uid, params)
        except Exception as exc:
            repo.change(uid, status="failed", error=str(exc), finished_at=now())
            raise
        return uid

    def execute(self, uid, params):
        app = self.app
        repo = app.extensions["jurassic_repo"]
        run_dir = Path(app.config["INSTANCE_DIR"]) / "runs" / uid
        root = Path(app.config["JURASSIC_ROOT"])
        start = time.monotonic()
        history = []
        logs_dir = run_dir / "logs"
        logs_dir.mkdir(exist_ok=True)
        def step(label, argv):
            repo.change(uid, status="running", stage=label)
            target = logs_dir / f"{label}.log"
            env = os.environ.copy()
            env["OMP_NUM_THREADS"] = str(app.config["OMP_THREADS"])
            dt = time.monotonic()
            process = subprocess.Popen(argv, cwd=run_dir, env=env, stdout=subprocess.PIPE,
                                       stderr=subprocess.STDOUT, start_new_session=True, text=True, errors="replace")
            try:
                out, _ = process.communicate(timeout=app.config["RUN_TIMEOUT"])
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                out, _ = process.communicate()
                target.write_text(out[-200000:])
                raise RuntimeError(f"{label} exceeded the execution timeout")
            target.write_text(out[-200000:])
            row = {"stage": label, "argv": argv, "seconds": round(time.monotonic()-dt, 3), "exit_code": process.returncode}
            history.append(row)
            (run_dir / "commands.json").write_text(json.dumps(history, indent=2))
            if process.returncode != 0: raise RuntimeError(f"{label} returned exit code {process.returncode}; see log")
            if any(p.is_file() and p.stat().st_size > app.config["MAX_RUN_BYTES"] for p in run_dir.rglob("*")):
                raise RuntimeError("Run exceeded maximum allowed output size")

        try:
            ds = datasets(app).get(params["dataset"])
            if not ds: raise RuntimeError("Selected dataset is no longer available")
            ctl = make_control(root, params, ds)
            (run_dir / "config.ctl").write_text(ctl)
            if not params.get("uploaded_atmosphere"):
                # Uses default JURASSIC climatology for known gases.
                step("climatology", [str(root / "src" / "climatology"), "config.ctl", "atm.tab"])
            step("geometry", [str(root / "src" / params["geometry"]), "config.ctl", "obs.tab"])
            step("formod", [str(root / "src" / "formod"), "config.ctl", "obs.tab", "atm.tab", "rad.tab"])
            info = summary(run_dir / "rad.tab", params["channels"])
            if info["rows"] > 500: raise RuntimeError("Output observation count exceeds GUI safety limit")
            if params["kernel"]:
                step("kernel", [str(root / "src" / "kernel"), "config.ctl", "obs.tab", "atm.tab", "kernel.tab"])
            repo.change(uid, stage="Exporting plots")
            (run_dir / "plots").mkdir(exist_ok=True)
            for metric in ("radiance", "transmittance", "atmosphere", "kernels"):
                if metric == "kernels" and not params["kernel"]:
                    continue
                try:
                    if metric in ("radiance", "transmittance"):
                        for ext in ("png", "svg"):
                            render_chart(run_dir / "rad.tab", params, metric,
                                         run_dir / "plots" / f"{metric}.{ext}")
                    elif metric == "atmosphere":
                        render_atmosphere(run_dir / "atm.tab", params["gases"],
                                          run_dir / "plots" / "atmosphere.png")
                    else:
                        render_kernel(run_dir / "kernel.tab", run_dir / "plots" / "kernels.png")
                except Exception as exc:
                    with (logs_dir / "plotting.log").open("a") as report:
                        report.write(f"{metric}: {exc}\n")
            warnings = sum((logs_dir / (stage+".log")).read_text(errors="replace").count("Warning (") for stage in ("formod", "kernel") if (logs_dir / (stage+".log")).exists())
            (run_dir / "metadata.json").write_text(json.dumps({"id": uid, "parameters":params, "results":info, "commands":history, "warnings": warnings},indent=2))
            repo.change(uid, status="completed", stage="Results ready", finished_at=now(),
                        duration=round(time.monotonic()-start,2), warnings=warnings)
        except Exception as exc:
            (logs_dir / "web_runner.log").write_text(traceback.format_exc()[-20000:])
            repo.change(uid, status="failed", stage="Execution failed", error=str(exc),
                        finished_at=now(), duration=round(time.monotonic()-start, 2))
