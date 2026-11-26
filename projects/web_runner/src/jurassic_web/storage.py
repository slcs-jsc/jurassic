"""SQLite metadata. Scientific files stay in private per-run folders."""
from __future__ import annotations
import json
import sqlite3
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path


def now():
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


class Repository:
    def __init__(self, filename):
        self.filename = Path(filename)

    @contextmanager
    def connect(self):
        db = sqlite3.connect(self.filename, timeout=20)
        db.row_factory = sqlite3.Row
        db.execute("PRAGMA journal_mode=WAL")
        try:
            yield db
            db.commit()
        except Exception:
            db.rollback()
            raise
        finally:
            db.close()

    def initialize(self):
        with self.connect() as db:
            db.execute("""CREATE TABLE IF NOT EXISTS runs (
                id TEXT PRIMARY KEY, created_at TEXT NOT NULL, finished_at TEXT,
                status TEXT NOT NULL, stage TEXT NOT NULL, geometry TEXT NOT NULL,
                dataset TEXT NOT NULL, parameters TEXT NOT NULL, error TEXT DEFAULT '',
                duration REAL DEFAULT 0, warnings INTEGER DEFAULT 0)""")

    def interrupt_stale_runs(self):
        with self.connect() as db:
            db.execute("UPDATE runs SET status='interrupted', stage='Stopped after service restart', finished_at=?, error='Worker interrupted during service restart' WHERE status IN ('queued','running')", (now(),))

    def create(self, id, parameters):
        with self.connect() as db:
            db.execute("INSERT INTO runs (id, created_at, status, stage, geometry, dataset, parameters) VALUES (?,?,?,?,?,?,?)",
                       (id, now(), "queued", "Waiting for worker", parameters["geometry"], parameters["dataset"], json.dumps(parameters)))

    def change(self, id, **fields):
        allowed = {"status", "stage", "error", "finished_at", "duration", "warnings"}
        update = {k:v for k,v in fields.items() if k in allowed}
        if not update:
            return
        keys = list(update)
        with self.connect() as db:
            db.execute("UPDATE runs SET " + ",".join(f"{k}=?" for k in keys) + " WHERE id=?", [update[k] for k in keys] + [id])

    @staticmethod
    def record(row):
        if row is None: return None
        d = dict(row)
        d["parameters"] = json.loads(d["parameters"])
        return d

    def get(self, id):
        with self.connect() as db:
            row = db.execute("SELECT * FROM runs WHERE id=?", (id,)).fetchone()
        return self.record(row)

    def recent(self, limit=40):
        with self.connect() as db:
            rows = db.execute("SELECT * FROM runs ORDER BY created_at DESC LIMIT ?", (limit,)).fetchall()
        return [self.record(row) for row in rows]

    def active_count(self):
        with self.connect() as db:
            return db.execute("SELECT COUNT(*) FROM runs WHERE status IN ('queued', 'running')").fetchone()[0]
