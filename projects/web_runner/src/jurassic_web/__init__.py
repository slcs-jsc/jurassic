"""Flask application factory. The web app never modifies upstream JURASSIC sources."""
from __future__ import annotations

import os
from pathlib import Path
from flask import Flask, abort, jsonify, request, session
from .storage import Repository
from .jobs import RunQueue
from .routes import bp
from .settings import SETTINGS


def create_app(test_config=None):
    app = Flask(__name__, template_folder="templates", static_folder="static")
    basedir = Path(__file__).resolve().parents[2]
    jurassic_root = Path(os.environ.get("JURASSIC_ROOT", str(basedir.parents[1]))).resolve()
    instance = Path(os.environ.get("JURASSIC_WEB_INSTANCE", str(basedir / "instance"))).resolve()
    app.config.update(
        JURASSIC_ROOT=jurassic_root,
        INSTANCE_DIR=instance,
        MAX_CONTENT_LENGTH=190 * 1024 * 1024,
        RUN_TIMEOUT=int(os.environ.get("JURASSIC_WEB_TIMEOUT", "180")),
        OMP_THREADS=max(1, min(8, int(os.environ.get("JURASSIC_WEB_OMP_THREADS", "2")))),
        MAX_ACTIVE_RUNS=6,
        MAX_RUN_BYTES=SETTINGS["max_run_bytes"],
        SESSION_COOKIE_HTTPONLY=True,
        SESSION_COOKIE_SAMESITE="Lax",
        SEND_FILE_MAX_AGE_DEFAULT=86400,
    )
    if test_config:
        app.config.update(test_config)
    instance = Path(app.config["INSTANCE_DIR"])
    instance.mkdir(parents=True, exist_ok=True)
    secretfile = instance / "flask_session_key"
    if not secretfile.exists():
        import secrets
        secretfile.write_bytes(secrets.token_bytes(32))
        secretfile.chmod(0o600)
    app.secret_key = secretfile.read_bytes()

    repo = Repository(instance / "runs.sqlite3")
    repo.initialize()
    repo.interrupt_stale_runs()
    app.extensions["jurassic_repo"] = repo
    app.extensions["jurassic_queue"] = RunQueue(app)

    @app.before_request
    def verify_csrf():
        import secrets
        if "csrf" not in session:
            session["csrf"] = secrets.token_urlsafe(32)
        if request.method in {"POST", "PUT", "PATCH", "DELETE"}:
            submitted = request.headers.get("X-CSRF-Token", "") or request.form.get("csrf", "")
            if not secrets.compare_digest(submitted, session["csrf"]):
                abort(403, "Invalid CSRF token. Reload the page and try again.")

    @app.context_processor
    def context():
        return {"csrf_token": session.get("csrf", ""), "app_name": "JURASSIC Web Runner"}

    @app.after_request
    def security_headers(response):
        response.headers["X-Content-Type-Options"] = "nosniff"
        response.headers["X-Frame-Options"] = "DENY"
        response.headers["Referrer-Policy"] = "same-origin"
        response.headers["Cache-Control"] = "no-store" if request.path.startswith("/api/") else "no-cache"
        response.headers["Content-Security-Policy"] = ("default-src 'self'; script-src 'self'; style-src 'self' 'unsafe-inline'; "
            "img-src 'self' data: blob:; connect-src 'self'; font-src 'self'; object-src 'none'; base-uri 'self'; frame-ancestors 'none'")
        return response

    @app.errorhandler(400)
    @app.errorhandler(403)
    @app.errorhandler(404)
    @app.errorhandler(413)
    def errors(exc):
        if request.path.startswith("/api/"):
            return jsonify(error=str(getattr(exc, "description", exc))), exc.code
        from html import escape
        return f"<h1>Request failed ({exc.code})</h1><p>{escape(str(exc.description))}</p>", exc.code

    app.register_blueprint(bp)
    return app
