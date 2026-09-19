#!/usr/bin/env python3
"""Run the complete single-core JURASSIC CGA validation spectrum calculation."""

import sys
from run_ega import add_common_arguments, execute, make_parser


def main():
    parser = make_parser(__doc__)
    add_common_arguments(parser)
    if len(sys.argv) == 1:
        parser.print_help()
        return
    args = parser.parse_args()
    try:
        execute("cga", "test_cga", args)
    except RuntimeError as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    main()
