#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

mapfile -t SESSIONS < <(tmux ls 2>/dev/null | awk -F: '/bench_.*(_num_features_extra|_all_features)/{print $1}')
for session in "${SESSIONS[@]}"; do
  tmux kill-session -t "$session" 2>/dev/null || true
done

NUM_FEATURES_SUPPLEMENT_START_OFFSET=0 \
NUM_FEATURES_SUPPLEMENT_MAX_JOBS=4 \
bash scripts/start_num_features_supplement_tmux.sh
