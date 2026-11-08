#!/usr/bin/env python3
"""Create synthetic truth from the unperturbed climatological retrieval prior."""

import math
import sys
from pathlib import Path

# Peak anomalies and their vertical locations; the common width keeps each
# perturbation resolved across several 2 km atmospheric grid levels.
TEMP_AMPLITUDE_K = 5
O3_AMPLITUDE_FRACTION = 0.25
HNO3_AMPLITUDE_FRACTION = -0.25
TEMP_CENTER_KM = 25
O3_CENTER_KM = 28
HNO3_CENTER_KM = 23
GAUSSIAN_WIDTH_KM = 5


def gaussian(z, center):
    """Return a smooth unit-height vertical anomaly at altitude z [km]."""
    return math.exp(-0.5 * ((z - center) / GAUSSIAN_WIDTH_KM) ** 2)


def main(source, target):
    # The input is the unperturbed a priori. Preserve every field except the
    # three retrieved quantities (pressure and background gases stay fixed).
    # JURASSIC ASCII columns: z=2, T=6, O3=9, HNO3=10 (1-based).
    rows = []
    for line in Path(source).read_text().splitlines():
        if not line or line.startswith("#"):
            rows.append(line)
            continue
        fields = line.split()
        z = float(fields[1])
        # Add temperature in K; apply O3/HNO3 as fractional VMR changes.
        # Smooth Gaussians make vertically resolved anomalies, not isolated levels.
        fields[5] = f"{float(fields[5]) + TEMP_AMPLITUDE_K * gaussian(z, TEMP_CENTER_KM):.10g}"
        fields[8] = f"{float(fields[8]) * (1 + O3_AMPLITUDE_FRACTION * gaussian(z, O3_CENTER_KM)):.10g}"
        fields[9] = f"{float(fields[9]) * (1 + HNO3_AMPLITUDE_FRACTION * gaussian(z, HNO3_CENTER_KM)):.10g}"
        rows.append(" ".join(fields))
    Path(target).write_text("\n".join(rows) + "\n")


if __name__ == "__main__":
    main(*sys.argv[1:])
