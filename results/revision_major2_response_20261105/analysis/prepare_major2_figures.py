from pathlib import Path
root=Path("/home/xd/y/Benchmark");rev=root/"results/revision_major2_response_20261105"
p=(root/"scripts/generate_fig4a_spatial_benchmark.R").read_text()
old='''    selected
}

methods_meta <- pick_representative_methods''';new='''    extra_random <- available |>
        filter(.data$MethodBase == "random", .data$Method == "random-N2000") |>
        mutate(IsBaseline = TRUE) |>
        select(.data$MethodBase, .data$Method, .data$IsBaseline)
    bind_rows(selected, extra_random)
}

methods_meta <- pick_representative_methods'''
assert p.count(old)==1;p=p.replace(old,new,1);(rev/"analysis/plot_fig4a_add_random2000.R").write_text(p)
p=(root/"scripts/generate_fig6_biological_context.R").read_text()
replacements=[('''full_scores_path <- file.path(
    "results", "spatial_svg_rebuild_v1", "frozen_scores",
    "dataset_representative_ranks.tsv"
)''','full_scores_path <- "results/revision_major2_response_20261105/fig5_inputs/full_ranks.tsv"'),('''subset_scores_path <- file.path(
    "results", "biological_subsets_rebuild_v1", "frozen_scores",
    "dataset_representative_ranks.tsv"
)''','subset_scores_path <- "results/revision_major2_response_20261105/fig5_inputs/subset_ranks.tsv"'),('''representatives_path <- file.path(
    "results", "biological_subsets_rebuild_v1", "frozen_scores",
    "representative_settings.tsv"
)''','representatives_path <- "results/revision_major2_response_20261105/fig5_inputs/representative_settings.tsv"'),('        base <- display_method(method)','        base <- if (method == "random_2000") "Random" else display_method(method)'),('''    method <- representatives$fs_method[[i]]
    n_value <- if''','''    method <- representatives$fs_method[[i]]
    if (method == "random_2000") next
    n_value <- if'''),('Expected = c(34, 34 * 3 * 3, 34 * 3, NA, NA, NA, 0, 0, 0)','Expected = c(35, 35 * 3 * 3, 34 * 3, NA, NA, NA, 0, 0, 0)')]
for old,new in replacements:
 assert p.count(old)==1,(old,p.count(old));p=p.replace(old,new,1)
(rev/"analysis/plot_fig5_add_random2000.R").write_text(p)
print("Prepared figure scripts with original themes and dimensions")
