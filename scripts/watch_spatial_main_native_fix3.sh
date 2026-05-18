#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

LOG_DIR="results/spatial_main_native_seed0_fix3/monitor"
LOG_FILE="$LOG_DIR/watch.log"
mkdir -p "$LOG_DIR"

status_cmd="./scripts/tmux_spatial_main_native_status.sh"

while true; do
  {
    echo "===== $(date '+%F %T %Z') ====="
    echo "[status]"
    "$status_cmd" || true
    echo
    echo "[latest-errors]"
    python - <<'PY'
from pathlib import Path

roots = sorted(Path("results/spatial_main_native_seed0_fix3").glob("*/logs/benchmark.log"))
for log_path in roots:
    lines = log_path.read_text(errors="ignore").splitlines() if log_path.exists() else []
    tail = lines[-80:]
    errors = [line for line in tail if "ERROR" in line or "Traceback" in line or "ValueError:" in line or "RuntimeError:" in line]
    if errors:
        print(f"--- {log_path}")
        for line in errors[-12:]:
            print(line)
PY
    echo
  } >> "$LOG_FILE" 2>&1
  sleep 1200
done
