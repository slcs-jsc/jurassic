"""Local development entrypoint. Use Waitress/appropriate reverse proxy for deployments."""
from __future__ import annotations
import argparse
import os
from . import create_app


def main():
    parser = argparse.ArgumentParser(description="Start JURASSIC Web Runner")
    parser.add_argument("--host", default=os.environ.get("JURASSIC_WEB_HOST", "127.0.0.1"))
    parser.add_argument("--port", type=int, default=int(os.environ.get("JURASSIC_WEB_PORT", "5000")))
    parser.add_argument("--serve", choices=["waitress", "flask"], default="waitress")
    args = parser.parse_args()
    app = create_app()
    if args.serve == "waitress":
        from waitress import serve
        serve(app, host=args.host, port=args.port, threads=4)
    else:
        app.run(host=args.host, port=args.port, debug=False, use_reloader=False)


if __name__ == "__main__":
    main()
