#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
OUT="results/resource_profile_v1"
mkdir -p "$OUT/json" "$OUT/time"

methods=(
    Brennecke TFs all_features anticor dubstepr hotspot morans_i nbumi nnsvg osca random
    scPNMF scanpy_cell_ranger scanpy_cell_ranger_batch scanpy_pearson scanpy_pearson_batch
    scanpy_seurat scanpy_seurat_batch scanpy_seurat_v3 scanpy_seurat_v3_batch scry
    scsegindex seurat_disp seurat_mvp seurat_sct seurat_vst singleCellHaystack somde
    sparkx spatialde statistic_mean statistic_variance triku wilcoxon
)

run_method() {
    local method="$1"
    if [[ -s "$OUT/json/${method}.json" ]]; then
        return
    fi
    /usr/bin/time -v -o "$OUT/time/${method}.txt" \
        .conda-env/bin/python scripts/profile_feature_selector.py \
        --method "$method" --output "$OUT/json/${method}.json" \
        >"$OUT/${method}.log" 2>&1
}

export -f run_method
export OUT ROOT
printf '%s\n' "${methods[@]}" | xargs -P 2 -I{} bash -c 'run_method "$@"' _ {}
