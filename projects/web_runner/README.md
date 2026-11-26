# JURASSIC Web Runner

A locally runnable Flask GUI for the upstream [JURASSIC](https://github.com/slcs-jsc/jurassic) scientific radiative-transfer executables. UI uses soft-neumorphism with **light mode default**, optional dark mode, responsive layouts and real interactive Plotly charts.

## Location and safety

Expected location: `~/jurassic/src/web-runner/`. The web runner is an independent Python package. It does **not** change JURASSIC's C sources, example `.ctl` files, or reference results. Runtime files and SQLite data are stored in `web-runner/instance/` (Git-ignored).

## Installation

See accompanying `install_jurassic_web_runner.sh`. It safely extracts the ZIP under `~/jurassic/src/web-runner/`, creates `.venv`, and uses pip to install this project and developer tools.

Start it:

```bash
cd ~/jurassic/src/web-runner
./start.sh
# Open http://127.0.0.1:5000
```

If your Python installation lacks venv support:

```bash
sudo apt install python3-venv python3-pip
```

Dependencies require network access to the configured package index on first install. The Plotly JavaScript asset is bundled locally and does not need a CDN. The matplotlib backend is headless (`Agg`).

To validate Python code:

```bash
./check.sh
```

## User interface

- **Simulation**: geometry (nadir, limb, zenith), dataset, available gas and channel selection, output type, ray-tracing parameters, optional kernels and custom atmosphere.
- **Configuration**: preview, download generated `.ctl`, import a `.ctl` with an explicit supported-key allowlist. Unknown control directives are ignored, and uploaded file paths are never executed.
- **Dataset Manager**: inspect lookup-table availability and register compatible extra gas/channel data as a ZIP of `.tab`/`.filt` files. This registers existing scientific data; it **does not** create missing spectroscopy.
- **Results**: live status, interactive charts, transmittance, atmospheric profile, optional sampled kernel heatmap, data table, logs, and configuration.
- **History & Comparison**: browse runs, open historical results, overlay 2–4 runs with matching geometry axes.
- **Downloads**: `rad.tab`, `atm.tab`, `obs.tab`, `config.ctl`, optional `kernel.tab`, publication-quality PNG/SVG, CSV, plots ZIP, full results ZIP.

No user accounts are present as requested. A CSRF session token protects state-changing requests. The app uses a built-in single-worker local queue and a SQLite database. Run statuses are persisted and interrupted jobs are identified after restart. The in-process queue is **development only** and must be replaced with a persistent external worker for reliable multi-process VM deployments.

## Dataset names

For example `boxcar_792.0000_CO2.tab` and `boxcar_792.0000.filt`. A ZIP dataset bundle must contain only flat files matching that pattern with a single family prefix. Each `*.tab` should contain compatible numerical spectral lookup-table data. Header-only tables are identified and warned about. Scientific validity **cannot** be established merely by parsing file names and numeric rows.

All selected channels require a filter and at least one numerical table; partial gas/channel support is shown as a warning. Custom gases require a matching uploaded atmospheric `.tab` file containing all selected gas names in its column headers. In addition, this version limits configured channels/gases (16 each), maximum returned observations (500), worker queue size (6), and execution duration (180 seconds per stage by default).

## Scientific scope

This software delegates **all radiative-transfer calculations** to upstream `climatology`, geometry executables, `formod`, and optional `kernel`. It does not fix or reinterpret JURASSIC model warnings. `make check` regression testing and example reference checks are distinct from successful simulation completion. In particular, the earlier Limb and Zenith reference discrepancies remain scientific diagnostics, not web GUI bugs.

## Project layout

```
src/jurassic_web/
    __init__.py         Flask factory and security headers
    routes.py           UI and API endpoints
    configuration.py    Constrained `.ctl` generation/import
    catalog.py          Lookup-table discovery and dataset registry
    jobs.py             Isolated command runner and status tracking
    storage.py          SQLite run metadata
    parsers.py          ASCII radiance, atmosphere and kernel readers
    plotting.py         Static Matplotlib exports
    templates/          Semantic, accessible HTML
    static/css/         Light/dark neumorphic design tokens and layout
    static/js/          Interactive client logic
    static/vendor/      Bundled Plotly JS
instance/               Runtime files, not committed
    runs/<id>/           Config, inputs, outputs, logs, plots, metadata
    datasets/<id>/       Registered user-supplied spectroscopy packages
tests/                  Unit and HTTP integration tests
```

## Before public VM deployment

This is a feature-rich **local working version**, not a production-hardened public multiuser service. Before making the VM accessible to everyone, add **reverse-proxy HTTPS, network-level or application request rate limits, disk quotas and cleanup, strong OS process isolation (container/unprivileged worker), a persistent queue (RQ/Celery), file antivirus/content scanning as appropriate, observability, and backup policies**. Running scientific executables on uploaded inputs is risky in a public anonymous setting. Do not simply launch the development service on `0.0.0.0` for public access without those controls.

## Configuration and environment

- `JURASSIC_ROOT`: path to upstream checkout (normally derived automatically)
- `JURASSIC_WEB_INSTANCE`: private state directory (default: `web-runner/instance`)
- `JURASSIC_WEB_HOST`: address to bind (default `127.0.0.1`)
- `JURASSIC_WEB_PORT`: listening port (default `5000`)
- `JURASSIC_WEB_TIMEOUT`: per-stage time limit in seconds (default `180`)
- `JURASSIC_WEB_OMP_THREADS`: OpenMP threads per simulation (default `2`, capped at `8`)

## Known limitations / next phase

- Control import intentionally supports a **subset** of the upstream keys; custom `.ctl` files are not executed verbatim.
- Uploaded spectral tables are checked for basic format and availability, not scientifically validated.
- Atmospheric uploads require JURASSIC-compatible columns; the app does not synthesize new gas profiles.
- In-process jobs do not survive application restart, and concurrent public VM hosting requires an external worker.
- ASCII `rad.tab` parser expects the standard `10 + 2×ND` layout; the initial UI targets example-compatible scenarios.
- Kernel heatmap is sampled for responsive rendering; the full `kernel.tab` remains downloadable.
