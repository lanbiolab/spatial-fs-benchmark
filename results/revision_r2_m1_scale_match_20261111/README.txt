Paired E9.5 aggregation sensitivity (Supplemental Figures S26-S28)

These scripts retain the executed analysis settings. Adapt ROOT=/home/xd/y/Benchmark
and dataset paths to the extracted benchmark installation; install the core and
method-specific environments described in the archive README. Run from benchmark/.

1. scripts/prepare.py constructs the two arms from the preprocessed E9.5 input and
   generates the configs and jobs.json. The two arms share the bead population,
   gene universe and per-gene totals. Nominal dimensions use the conditional 0.65
   source-coordinate convention; this is not an independently calibrated assay.
2. scripts/run_job.py CONFIG runs one job using the existing benchmark wrappers.
   Iterate over the configurations in jobs.json, adapting device/thread resources
   for the local machine. The archived dispatch_jobs.py documents the original
   phased scheduler; it expects completion gates from the historical pilot run.
   The original paired comparison records the SPARK-X pooled failure explicitly;
   missing runs must not be substituted with zero or worst scores.
3. scripts/analyze.py and verify.py generate and validate the per-seed metrics,
   actual-count checks, endpoint ranks and rank correlations after run completion.
   documented_method_failures.json records the diagnosed implementation exception.
4. scripts/export_geometry.py and plot.R generate geometry and complete rank
   panels. The final summary layout is ../revision_r2_m1_response_20261111/scripts/
   plot_summary.R; copy the generated analysis inputs to that directory's analysis/
   folder or adapt its input root. Fonts and the original shared plotting helpers
   are described in the archive README. The figures are reused in S26-S28 without
   changing data or graphical mappings.

No primary expression matrices or full derived-result tables are bundled here.
Public input accessions are in Methods; detailed analysis records are available
on request. This extension does not change the five-dataset primary ranking.
