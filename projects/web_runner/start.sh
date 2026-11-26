#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
if [[ ! -x "$ROOT/.venv/bin/jurassic-web" ]]; then
  echo "Virtual environment missing. Run the installer first, or: python3 -m venv '$ROOT/.venv' && '$ROOT/.venv/bin/pip' install -e '$ROOT'" >&2
  exit 1
fi
export JURASSIC_ROOT="${JURASSIC_ROOT:-$(cd "$ROOT/../.." && pwd)}"
exec "$ROOT/.venv/bin/jurassic-web" "$@"
