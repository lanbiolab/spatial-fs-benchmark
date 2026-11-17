set -e
cd /home/xd/y/Benchmark
export PYTHONPATH="$PWD/src"
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 NUMEXPR_NUM_THREADS=2
.conda-env/bin/python results/revision_r2_major4_20261107/code/run_retention.py > results/revision_r2_major4_20261107/run.log 2>&1
tail -n 12 results/revision_r2_major4_20261107/run.log
