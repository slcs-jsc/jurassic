#!/usr/bin/env python3
"""Generate the complete channel-matched RFM reference spectra."""

import os
import sys
from pathlib import Path
from run_ega import add_common_arguments, execute, make_parser


def main():
    parser = make_parser(__doc__)
    add_common_arguments(parser)
    parser.add_argument("--rfm-bin", type=Path,
                        default=os.environ.get("RFM_BIN", Path.home() / "wrk/rfm/v521_timings/rfm"),
                        help="instrumented external RFM executable (environment: RFM_BIN)")
    parser.add_argument("--rfm-hit", type=Path,
                        default=os.environ.get("RFM_HIT", Path.home() / "wrk/rfm/hitbin20/hitran2020_mir.bin"),
                        help="external HITRAN binary line file (environment: RFM_HIT)")
    parser.add_argument("--rfm-xsc-dir", type=Path,
                        default=os.environ.get("RFM_XSC_DIR", Path.home() / "wrk/rfm/xsc20"),
                        help="external RFM cross-section directory (environment: RFM_XSC_DIR)")
    if len(sys.argv) == 1:
        parser.print_help()
        return
    args = parser.parse_args()
    try:
        execute("rfm", "rfm_reference", args,
                ("--rfm-bin", str(args.rfm_bin), "--rfm-hit", str(args.rfm_hit),
                 "--rfm-xsc-dir", str(args.rfm_xsc_dir), "--require-rfm-timers"))
    except RuntimeError as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    main()
