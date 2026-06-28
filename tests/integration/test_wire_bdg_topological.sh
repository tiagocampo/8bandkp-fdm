#!/bin/bash
# U8 regression: wire BdG open->close->reopen + no auto-window fallback.
# Args: <topologicalAnalysis_exe> <config_file>
set -euo pipefail
EXE="$(realpath "$1")"
CONFIG="$(realpath "$2")"
HERE="$(cd "$(dirname "$0")" && pwd)"
python3 "$HERE/verify_wire_bdg_topological.py" "$EXE" "$CONFIG"
