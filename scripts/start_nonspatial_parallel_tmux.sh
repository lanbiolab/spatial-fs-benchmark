#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

launch_worker() {
    local gpu="$1"
    local worker="$2"
    shift 2
    local session="spatial_nonspatial_${worker}"
    tmux new-session -d -s "$session" \
        "cd '$ROOT' && bash scripts/run_nonspatial_parallel_worker.sh '$gpu' '$worker' $*"
}

launch_worker 0 gpu0 \
    configs/rebuild_v1/non_spatial_downstream/canonical/stomics_0212.yaml \
    configs/rebuild_v1/non_spatial_downstream/canonical/e8p5_embryo.yaml \
    configs/rebuild_v1/non_spatial_downstream/feature_number/stomics_0212.yaml \
    configs/rebuild_v1/non_spatial_downstream/feature_number/e8p5_embryo.yaml

launch_worker 1 gpu1 \
    configs/rebuild_v1/non_spatial_downstream/canonical/stomics_0218.yaml \
    configs/rebuild_v1/non_spatial_downstream/canonical/e9p5_embryo.yaml \
    configs/rebuild_v1/non_spatial_downstream/feature_number/stomics_0218.yaml \
    configs/rebuild_v1/non_spatial_downstream/feature_number/e9p5_embryo.yaml

launch_worker 2 gpu2 \
    configs/rebuild_v1/non_spatial_downstream/canonical/dlpfc.yaml \
    configs/rebuild_v1/non_spatial_downstream/canonical/stomics_0224.yaml \
    configs/rebuild_v1/non_spatial_downstream/feature_number/mouse_brain_serial_sections.yaml \
    configs/rebuild_v1/non_spatial_downstream/feature_number/dlpfc.yaml \
    configs/rebuild_v1/non_spatial_downstream/feature_number/stomics_0224.yaml
