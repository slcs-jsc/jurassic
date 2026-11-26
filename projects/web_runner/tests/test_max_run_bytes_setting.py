"""Regression test for MAX_RUN_BYTES Flask configuration."""
from jurassic_web.settings import SETTINGS

def test_max_run_bytes_defined(app):
    assert app.config["MAX_RUN_BYTES"] == SETTINGS["max_run_bytes"]
