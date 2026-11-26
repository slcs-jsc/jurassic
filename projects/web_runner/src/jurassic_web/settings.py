"""UI and execution constraints; not JURASSIC scientific defaults."""
SETTINGS = {
    "geometries": ("nadir", "limb", "zenith"),
    "default_family": {"nadir": "airs", "limb": "boxcar", "zenith": "boxcar"},
    "safe_gases": ("CO2", "H2O", "O3", "F11", "CCl4"),
    "max_channels": 16,
    "max_gases": 16,
    "max_obs": 500,
    "max_run_bytes": 200 * 1024 * 1024,
}
