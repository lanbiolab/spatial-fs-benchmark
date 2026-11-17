#!/usr/bin/env bash
set -Eeuo pipefail
cd /home/xd/y/Benchmark
revision=results/revision_major2_random2000_20261105
mkdir -p "$revision/logs"
exec 9>/tmp/spatial_fs_benchmark_gpu_3.lock
if ! flock -n 9; then echo "FAILED GPU 3 lock is held"; touch "$revision/FAILED"; exit 1; fi
if nvidia-smi -i 3 --query-compute-apps=pid --format=csv,noheader,nounits | grep -q '[0-9]'; then echo "FAILED GPU 3 already has another compute process"; touch "$revision/FAILED"; exit 1; fi
export CUDA_VISIBLE_DEVICES=3 OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 R_FUTURE_PLAN=sequential
run_one() {
 local name="$1"
 echo "[$(date -Is)] START $name" | tee -a "$revision/logs/status.log"
 if .conda-env/bin/spatial-fs-run --config "$revision/configs/$name.yaml" >"$revision/logs/$name.log" 2>&1 && .conda-env/bin/python "$revision/scripts/validate_random2000_run.py" "$name" >"$revision/logs/${name}_validation.log" 2>&1; then
  echo "[$(date -Is)] VALIDATED $name" | tee -a "$revision/logs/status.log"
 else
  code=$?; echo "[$(date -Is)] FAILED $name code=$code" | tee -a "$revision/logs/status.log"; return "$code"
 fi
}
worker_a() { for name in stomics_0224; do run_one "$name" || return; done; }
worker_b() { for name in mouse_brain_serial_sections dlpfc; do run_one "$name" || return; done; }
echo "[$(date -Is)] QUEUE RESUME GPU=3 workers=2 after verified legacy metadata absence" | tee -a "$revision/logs/status.log"
worker_a & a=$!
worker_b & b=$!
sa=0;sb=0;wait "$a" || sa=$?;wait "$b" || sb=$?
if (( sa != 0 || sb != 0 )); then echo "[$(date -Is)] QUEUE FAILED workerA=$sa workerB=$sb" | tee -a "$revision/logs/status.log";touch "$revision/FAILED";exit 1;fi
touch "$revision/DONE"
echo "[$(date -Is)] QUEUE DONE" | tee -a "$revision/logs/status.log"
