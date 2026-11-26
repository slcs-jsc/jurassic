"""Flask UI and JSON endpoints; all inputs validated server-side."""
from __future__ import annotations
import io
import json
import os
import re
import threading
import shutil
import uuid
import zipfile
from pathlib import Path
from flask import Blueprint, Response, abort, current_app as app, jsonify, render_template, request, send_file, send_from_directory
from werkzeug.utils import secure_filename
from .catalog import NAME, compatible, datasets, scan_family
from .configuration import import_form, make_control, normalize
from .parsers import atmospheric_chart, chart, rows_preview, summary, kernel_heatmap
from .settings import SETTINGS

bp = Blueprint("main", __name__)
VALID_RUN = re.compile(r"^[a-f0-9]{16}$")
FILE_NAMES = {"config.ctl", "rad.tab", "atm.tab", "obs.tab", "kernel.tab", "commands.json", "metadata.json"}
PLOT_METRICS = frozenset({"radiance", "transmittance", "atmosphere", "kernels"})
PLOT_LOCK = threading.RLock()


def repo(): return app.extensions["jurassic_repo"]


def output_dir(uid):
    if not VALID_RUN.fullmatch(uid): abort(404)
    if not repo().get(uid): abort(404)
    return Path(app.config["INSTANCE_DIR"]) / "runs" / uid


def issue(message, code=400): return jsonify({"error":str(message)}), code


@bp.route("/")
def index():
    return render_template("dashboard.html", initial_runs=repo().recent(8))


@bp.route("/api/system")
def system():
    root = Path(app.config["JURASSIC_ROOT"])
    required = ("climatology", "nadir", "limb", "zenith", "formod", "kernel")
    status = {name: (root / "src" / name).is_file() and os.access(root / "src" / name, os.X_OK) for name in required}
    revision = "unknown"
    try:
        import subprocess
        r = subprocess.run(["git", "-C", str(root), "rev-parse", "--short", "HEAD"],
                           capture_output=True, text=True, timeout=3)
        if r.returncode == 0: revision = r.stdout.strip()
    except (OSError, subprocess.TimeoutExpired): pass
    return jsonify({"ready": all(status.values()), "executables": status, "revision":revision,
                    "datasets":len(datasets(app)), "queued":repo().active_count()})


@bp.route("/api/catalog")
def catalog():
    items = datasets(app)
    for item in items.values(): item.pop("prefix", None)  # filesystem paths not needed by clients
    return jsonify(items)


@bp.route("/api/catalog/upload", methods=["POST"])
def upload_catalog():
    attachment = request.files.get("dataset")
    label = str(request.form.get("label", "")).strip()[:55]
    if not attachment or not attachment.filename.lower().endswith(".zip"): return issue("Upload a ZIP containing .tab and .filt files")
    if not (1 <= len(label) <= 55) or any(c in label for c in "<>\\/"):
        return issue("Provide a short dataset label")
    try:
        data = attachment.read(150 * 1024 * 1024 + 1)
        if len(data) > 150 * 1024 * 1024: return issue("Dataset ZIP exceeds 150 MB",413)
        archive = zipfile.ZipFile(io.BytesIO(data))
        entries = [e for e in archive.infolist() if not e.is_dir()]
        if not entries or len(entries) > 64: return issue("ZIP must contain 1–64 files")
        if sum(e.file_size for e in entries) > 220 * 1024 * 1024: return issue("Decompressed data exceeds 220 MB")
        families = set()
        for e in entries:
            # Prevent arbitrary paths, duplicate filenames, symlinks and unknown extensions.
            if e.filename != Path(e.filename).name or len(e.filename) > 100 or e.file_size > 22 * 1024 * 1024:
                return issue("ZIP contains an unsafe filename or oversized entry")
            if (e.external_attr >> 16) & 0o170000 == 0o120000: return issue("ZIP symlinks are not allowed")
            match = NAME.fullmatch(e.filename)
            if not match: return issue("Unsupported file name: " + e.filename)
            families.add(match["family"])
        if len(families) != 1: return issue("ZIP must contain one dataset family prefix")
        if len(set(e.filename for e in entries)) != len(entries): return issue("Duplicate filenames in ZIP")
        slug = "user-" + uuid.uuid4().hex[:12]
        directory = Path(app.config["INSTANCE_DIR"]) / "datasets" / slug
        directory.mkdir(parents=True, exist_ok=False)
        try:
            for e in entries:
                with archive.open(e) as src, (directory / e.filename).open("wb") as dest:
                    shutil.copyfileobj(src,dest,1024*1024)
            fam = families.pop()
            ds = scan_family(directory,fam)
            if not any(v["filter"] and any(g["data"] for g in v["gases"].values()) for v in ds["channels"].values()):
                raise ValueError("Dataset needs a filter and numeric gas table for at least one channel")
            (directory / "manifest.json").write_text(json.dumps({"family":fam,"label":label},indent=2))
        except Exception:
            shutil.rmtree(directory)
            raise
        return jsonify({"id":slug,"label":label,"channels":len(ds["channels"]), "message":"Dataset registered; check scientific compatibility before use."}),201
    except (zipfile.BadZipFile, RuntimeError, ValueError, OSError) as exc:
        return issue(str(exc))


def get_params(payload):
    if not isinstance(payload, dict): raise ValueError("Expected a JSON object")
    dataset_id = str(payload.get("dataset", ""))
    ds = datasets(app).get(dataset_id)
    if not ds: raise ValueError("Select an available spectroscopic dataset")
    params = normalize(payload, ds)
    return params, ds


@bp.route("/api/config/import",methods=["POST"])
def import_config():
    file = request.files.get("configuration")
    if not file or not file.filename.lower().endswith(".ctl"): return issue("Upload a .ctl configuration file")
    contents = file.read(100001)
    try:
        data = import_form(contents.decode("utf-8"))
        return jsonify({"values":data, "note":"Only supported UI parameters were imported; paths and unrecognized directives were ignored."})
    except (UnicodeDecodeError, ValueError) as exc: return issue(str(exc))


@bp.route("/api/config/preview",methods=["POST"])
def preview_config():
    try:
        params, ds = get_params(request.get_json(silent=True) or {})
        output = make_control(app.config["JURASSIC_ROOT"], params, ds)
        return jsonify({"text":output, "warnings":compatible(ds, params["channels"], params["gases"])})
    except ValueError as exc: return issue(str(exc))


@bp.route("/api/runs", methods=["GET"])
def history():
    return jsonify(repo().recent(100))


@bp.route("/api/runs", methods=["POST"])
def create_run():
    try:
        if request.mimetype == "multipart/form-data":
            payload = json.loads(request.form.get("payload", "{}"))
            file = request.files.get("atmosphere")
            if file and not file.filename.lower().endswith(".tab"):
                raise ValueError("Atmospheric upload must be a .tab file")
            atmosphere = file.read(3*1024*1024+1) if file else None
        else:
            payload = request.get_json(silent=True) or {}
            atmosphere = None
        params, ds = get_params(payload)
        if any(g not in SETTINGS["safe_gases"] for g in params["gases"]) and atmosphere is None:
            raise ValueError("Custom gases require an uploaded atmospheric profile containing matching gas columns")
        for exe in ("climatology",params["geometry"], "formod") + (("kernel",) if params["kernel"] else ()):
            binary = Path(app.config["JURASSIC_ROOT"]) / "src" / exe
            if not binary.is_file() or not os.access(binary,os.X_OK):
                raise ValueError(f"JURASSIC executable missing: {exe}. Build upstream with cd src && make")
        uid = app.extensions["jurassic_queue"].submit(params, atmosphere)
        return jsonify({"id": uid, "status":"queued", "warnings":compatible(ds,params["channels"],params["gases"])}), 202
    except (ValueError, TypeError, json.JSONDecodeError) as exc: return issue(str(exc))


@bp.route("/api/runs/<uid>")
def run_detail(uid):
    if not VALID_RUN.fullmatch(uid): abort(404)
    item = repo().get(uid)
    if not item: abort(404)
    if item["status"] == "completed":
        p = output_dir(uid) / "metadata.json"
        if p.is_file():
            item["results"] = json.loads(p.read_text()).get("results", {})
    return jsonify(item)


@bp.route("/api/runs/<uid>/chart/<metric>")
def run_chart(uid,metric):
    if metric not in ("radiance","transmittance","atmosphere"): abort(404)
    item = repo().get(uid)
    if not item: abort(404)
    if item["status"] != "completed": return issue("Run results not ready", 409)
    path = output_dir(uid)
    try:
        if metric == "atmosphere":
            return jsonify(atmospheric_chart(path / "atm.tab", item["parameters"]["gases"]))
        return jsonify(chart(path / "rad.tab", item["geometry"], item["parameters"]["channels"], metric))
    except ValueError as exc: return issue(str(exc),422)


@bp.route("/api/runs/<uid>/plot/<metric>.png")
def run_plot_png(uid, metric):
    """Show an existing PNG, or render a missing plot from stored scientific files.

    Completed historical runs are supported: old ZIPs contain radiance and
    transmittance PNGs but may not contain atmosphere/kernel PNGs yet.
    """
    if metric not in PLOT_METRICS: abort(404)
    theme = request.args.get("theme", "light")
    if theme not in ("light", "dark"):
        return issue("Unknown plot theme")
    run = repo().get(uid) if VALID_RUN.fullmatch(uid) else None
    if run is None: abort(404)
    if run["status"] != "completed": return issue("Run results not ready", 409)
    if metric == "kernels" and not run["parameters"].get("kernel"):
        return issue("Kernel calculation was not enabled for this run", 404)

    run_dir = output_dir(uid)
    plots = run_dir / "plots"
    filename = metric + ("_dark" if theme == "dark" else "") + ".png"
    image = plots / filename
    if not image.is_file():
        with PLOT_LOCK:
            if not image.is_file():
                from .plotting import render_chart, render_atmosphere, render_kernel
                plots.mkdir(parents=True, exist_ok=True)
                temporary = plots / (".rendering_" + filename)
                try:
                    if metric in ("radiance", "transmittance"):
                        render_chart(run_dir / "rad.tab", run["parameters"], metric, temporary, theme)
                    elif metric == "atmosphere":
                        render_atmosphere(run_dir / "atm.tab", run["parameters"]["gases"], temporary, theme)
                    else:
                        render_kernel(run_dir / "kernel.tab", temporary, theme)
                    temporary.replace(image)
                except (OSError, ValueError, RuntimeError, KeyError) as exc:
                    temporary.unlink(missing_ok=True)
                    return issue("Plot generation failed: " + str(exc), 422)
    return send_from_directory(plots, filename, mimetype="image/png",
                               as_attachment=request.args.get("download") == "1",
                               download_name=f"jurassic_{uid}_{filename}")


@bp.route("/api/runs/<uid>/table")
def run_table(uid):
    item = repo().get(uid)
    if not item: abort(404)
    if item["status"] != "completed": return issue("Run results not ready",409)
    try: return jsonify(rows_preview(output_dir(uid) / "rad.tab"))
    except ValueError as exc: return issue(str(exc),422)


@bp.route("/api/runs/<uid>/log")
def run_log(uid):
    output_dir(uid)
    step = request.args.get("step", "formod")
    if step not in ("climatology","geometry","formod","kernel","web_runner","plotting"):
        abort(404)
    file = output_dir(uid) / "logs" / (step + ".log")
    if not file.is_file(): return issue("Log not available",404)
    return jsonify({"step":step,"text":file.read_text(errors="replace")[-100000:]})


@bp.route("/api/runs/<uid>/download/<name>")
def download(uid,name):
    folder = output_dir(uid)
    if name == "rad.csv":
        import csv
        from .parsers import read_table
        data, headers = read_table(folder / "rad.tab")
        out = io.StringIO()
        writer = csv.writer(out)
        writer.writerow(headers or [f"Column {i+1}" for i in range(data.shape[1])])
        writer.writerows(data.tolist())
        return Response(out.getvalue(), mimetype="text/csv",
                        headers={"Content-Disposition": f'attachment; filename="jurassic_{uid}_rad.csv"'})
    if name in FILE_NAMES:
        path = folder / name
        if not path.is_file(): abort(404)
        return send_from_directory(folder,name,as_attachment=True)
    if name in {"radiance.png","transmittance.png","radiance.svg","transmittance.svg"}:
        return send_from_directory(folder / "plots", name, as_attachment=True)
    if name in ("all_results.zip", "plots.zip"):
        archive_path = folder / name
        with zipfile.ZipFile(archive_path,"w",compression=zipfile.ZIP_DEFLATED, compresslevel=6) as z:
            for file in folder.rglob("*"):
                if file.is_file() and file.suffix != ".zip":
                    if name == "plots.zip" and "plots" not in file.parts: continue
                    z.write(file,file.relative_to(folder))
        return send_file(archive_path,as_attachment=True,download_name=f"jurassic_{uid}_{name}")
    abort(404)


@bp.route("/api/compare", methods=["GET"])
def compare():
    ids = request.args.get("ids", "").split(",")
    metric = request.args.get("metric", "radiance")
    if len(ids) < 2 or len(ids) > 4 or metric not in ("radiance","transmittance"):
        return issue("Choose 2–4 runs and a supported metric")
    traces = []
    xlabel = ylabel = ""
    for uid in ids:
        item = repo().get(uid) if VALID_RUN.fullmatch(uid) else None
        if not item or item["status"] != "completed": return issue("All compared runs must be completed",409)
        data = chart(output_dir(uid) / "rad.tab",item["geometry"],item["parameters"]["channels"],metric)
        if traces and xlabel != data["xlabel"]: return issue("Comparison requires the same geometry axis")
        xlabel, ylabel = data["xlabel"],data["ylabel"]
        traces.extend({"name":f"{uid[:6]} · {t['name']}", "x":data["x"], "y":t["y"]} for t in data["traces"])
    return jsonify({"xlabel":xlabel,"ylabel":ylabel,"traces":traces})


@bp.route("/api/runs/<uid>/kernel")
def kernel_data(uid):
    item = repo().get(uid)
    if not item: abort(404)
    file = output_dir(uid) / "kernel.tab"
    if item["status"] != "completed" or not file.is_file(): return issue("Kernel output not available", 409)
    quantity = request.args.get("quantity", "TEMPERATURE")
    if quantity not in ("TEMPERATURE", "PRESSURE") and quantity not in item["parameters"]["gases"]:
        return issue("Unsupported sensitivity quantity")
    try: return jsonify(kernel_heatmap(file, quantity=quantity))
    except ValueError as exc: return issue(str(exc), 422)
