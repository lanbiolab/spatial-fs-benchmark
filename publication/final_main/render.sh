set -e
cd /home/xd/y/Benchmark
cp -a results/revision_arial_unified_20261109/inputs results/revision_arial_spacing_20261109/
FONTCONFIG_FILE=/home/xd/y/Benchmark/results/revision_arial_spacing_20261109/fonts.conf LD_PRELOAD=/home/xd/y/Benchmark/results/revision_arial_spacing_20261109/fractional_text.so .conda-env/bin/Rscript results/revision_arial_spacing_20261109/render.R > results/revision_arial_spacing_20261109/render.log 2>&1 || { tail -n 40 results/revision_arial_spacing_20261109/render.log; exit 1; }
cd results/revision_arial_spacing_20261109
zip -q plots.zip Figure*.pdf font_data_verification.txt render.log fractional_text.c
tail -n 8 render.log
