set -e
cd /home/xd/y/Benchmark
python3 - <<'PY'
from zipfile import ZipFile
ZipFile('results/revision_all_typography_20261109/inputs.zip').extractall('results/revision_all_typography_20261109/inputs')
PY
FONTCONFIG_FILE=/home/xd/y/Benchmark/results/revision_arial_spacing_20261109/fonts.conf LD_PRELOAD=/home/xd/y/Benchmark/results/revision_arial_spacing_20261109/fractional_text.so .conda-env/bin/Rscript results/revision_all_typography_20261109/render.R > results/revision_all_typography_20261109/render.log 2>&1 || { tail -n 50 results/revision_all_typography_20261109/render.log; exit 1; }
cd results/revision_all_typography_20261109
zip -q plots.zip Supplemental_Fig*.pdf font_data_verification.txt render.log
tail -n 5 render.log
