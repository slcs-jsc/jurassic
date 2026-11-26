# Design notes

## Boundary

Browser → Flask API → validated simulation specifications → isolated working folder → existing C executables → parsed scientific files → plots and downloads.

## Security

No shell interpolation. Fixed executable allowlist and isolated working directory. Programmatic generated control paths are sourced from the installed dataset catalog, not uploaded control-file directives. Uploaded `.ctl` parameters are allowlisted. File downloads are restricted to a fixed list and a validated run ID. Dataset bundles reject path traversal and symlinks. Session CSRF token is required on all POST requests. Per-stage timeout and output-size caps are applied.

## Scientifically relevant constraints

The existing lookup tables define the combinations the UI can offer. The app distinguishes header-only data and unregistered channel/gas files. Climatology is only used with known built-in emitters; custom gases require uploaded atmospheric profiles. The app records warnings and C process exit codes rather than suppressing them.

## VM readiness

See README deployment precautions. Most urgently: move runs out of Flask process into a supervised worker and provide application-level workload throttling and OS-level isolation before opening anonymous access from the internet.
