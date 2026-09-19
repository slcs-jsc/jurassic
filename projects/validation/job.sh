#!/usr/bin/env bash

set -euo pipefail

# Local defaults for the external validation data.  Exporting any of these
# variables before calling this script overrides the corresponding default.
export JURASSIC_TBL_DIR="${JURASSIC_TBL_DIR:-$HOME/wrk/jurassic/tab/tria_1cm/nc_1e-6}"
export RFM_BIN="${RFM_BIN:-$HOME/wrk/rfm/v521_timings/rfm}"
export RFM_HIT="${RFM_HIT:-$HOME/wrk/rfm/hitbin20/hitran2020_mir.bin}"
export RFM_XSC_DIR="${RFM_XSC_DIR:-$HOME/wrk/rfm/xsc20}"

# Each forward model remains single threaded.  Two independent 128-channel
# chunks can run concurrently on the reference notebook's performance cores.
export OMP_NUM_THREADS=1
export VALIDATION_JOBS="${VALIDATION_JOBS:-2}"
export VALIDATION_CPUSET="${VALIDATION_CPUSET:-0,2}"

root=$(cd "$(dirname "$0")/../.." && pwd)

[[ -d "$JURASSIC_TBL_DIR" ]] || {
  echo "Missing JURASSIC lookup-table directory: $JURASSIC_TBL_DIR" >&2
  exit 1
}
[[ -x "$RFM_BIN" ]] || {
  echo "Missing or non-executable RFM binary: $RFM_BIN" >&2
  exit 1
}
[[ -f "$RFM_HIT" ]] || {
  echo "Missing HITRAN file: $RFM_HIT" >&2
  exit 1
}
[[ -d "$RFM_XSC_DIR" ]] || {
  echo "Missing RFM cross-section directory: $RFM_XSC_DIR" >&2
  exit 1
}

# Remove incomplete chunk data and logs from earlier or interrupted runs.  The
# compact reference results remain in place until each new run has succeeded.
rm -rf "$root/projects/validation/work"

# Build with enough gas slots for the fixed 36-gas validation setup.
make -C "$root/src" clean
make -C "$root/src" DEFINES=-DNG=40 -j

# Regenerate compact reference spectra, approximation results, timings, and
# the complete report.  The Python runners retain the 128-channel chunking.
python3 "$root/projects/validation/run_rfm.py" --force
python3 "$root/projects/validation/run_ega.py" --force
python3 "$root/projects/validation/run_cga.py" --force
python3 "$root/projects/validation/analyze.py"

# A successful run has copied all compact results into their final locations.
rm -rf "$root/projects/validation/work"
