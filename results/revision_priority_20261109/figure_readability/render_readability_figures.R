#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(ggplot2)
    library(patchwork)
    library(ggforce)
    library(colorspace)
})

out <- 'results/revision_priority_20261109/figure_readability'
dir.create(out, recursive = TRUE, showWarnings = FALSE)
fig3 <- readRDS('results/revision_minor1_3_20261106/figures/figure_num_features_benchmark.rds')
fig4a <- readRDS('results/revision_major2_response_20261105/figures/fig4a/figure_fig4a_spatial_benchmark.rds')
fig5 <- readRDS('results/revision_major2_response_20261105/figures/fig5/figure_fig6_biological_context.rds')
fig6 <- readRDS('revision_r1_major1_five_20261105/figures/fig6/figure_fig6_spatial_integration.rds')

# Reuse the saved plot objects made by the original benchmark plotting scripts.
# Change the presentation themes and place each dense block on its own page.
readable_theme <- theme(
    axis.text.x = element_text(size = 7.8),
    axis.text.y = element_text(size = 7.8),
    axis.title = element_text(size = 8.2),
    strip.text = element_text(size = 8.2),
    plot.title = element_text(size = 8.6),
    plot.tag = element_text(size = 9.5, face = 'bold'),
    legend.text = element_text(size = 7.8),
    legend.title = element_text(size = 8.1),
    legend.key.width = grid::unit(13, 'pt'),
    legend.key.height = grid::unit(10, 'pt'),
    legend.spacing.x = grid::unit(3, 'pt')
)

save_figure <- function(name, plot, height) {
    # Letter paper with 1-inch margins gives a 6.5-inch column at native size.
    set.seed(20261109)
    ggsave(file.path(out, name), plot, width = 6.5, height = height,
           units = 'in', device = grDevices::cairo_pdf, family = 'Arial', bg = 'white')
}

assert_data_identity <- function(original, revised) {
    stopifnot(identical(original$data, revised$data),
              identical(original$mapping, revised$mapping),
              length(original$layers) == length(revised$layers))
    for (j in seq_along(original$layers)) {
        stopifnot(identical(original$layers[[j]]$layer_data(original$data),
                            revised$layers[[j]]$layer_data(revised$data)),
                  identical(original$layers[[j]]$mapping, revised$layers[[j]]$mapping))
    }
}

fig3_main <- fig3 & readable_theme
fig3_main[[1]] <- fig3_main[[1]] + theme(
    axis.title.y.left = element_text(vjust = 0.5, margin = margin(r = 5)),
    axis.text.x = element_blank(), strip.text = element_text(size = 8.2)
)
fig3_main[[2]] <- fig3_main[[2]] + theme(axis.text.x = element_blank(), strip.text = element_blank())
fig3_main[[3]] <- fig3_main[[3]] + theme(strip.text = element_blank())
for (i in 1:3) assert_data_identity(fig3[[i]], fig3_main[[i]])
save_figure('Figure3_readable.pdf', fig3_main, 7.2)

fig4a_companion <- fig4a & readable_theme
fig4a_companion[[2]] <- fig4a_companion[[2]] + theme(axis.text.y = element_blank())
rank_guide <- fig4a[[2]]$scales$get_scales('alpha')$clone()
stopifnot(identical(as.numeric(rank_guide$breaks), c(10, 20, 35)),
          identical(rank_guide$labels, c('10', '20', '30')))
rank_guide$labels <- as.character(rank_guide$breaks)
fig4a_companion[[2]] <- fig4a_companion[[2]] + rank_guide
fig4a_companion <- fig4a_companion + plot_annotation(tag_levels = NULL) & theme(
    legend.box = 'vertical', legend.title.position = 'top',
    legend.box.just = 'left', legend.position = 'bottom'
)
for (i in 1:2) assert_data_identity(fig4a[[i]], fig4a_companion[[i]])
save_figure('Supplemental_Figure_S10a.pdf', fig4a_companion, 6.6)

names6 <- c('b', 'c', 'd', 'e')
labels6 <- c('Mean scores', 'Score differences to scVI', 'Mean ranks', 'Rank differences to scVI')
fig6_block <- fig6[[1]]
stopifnot(length(fig6_block) == 4L)
for (i in seq_along(names6)) {
    p <- fig6_block[[i]] + readable_theme + labs(tag = NULL, title = labels6[[i]], x = NULL, y = NULL) + theme(
        axis.text.y = element_text(size = 7.8, face = 'plain', colour = 'black'),
        axis.ticks.y = element_line(), legend.position = 'bottom',
        legend.title.position = 'top', legend.box = 'horizontal',
        plot.margin = margin(5, 8, 5, 5), strip.text = element_text(size = 8.2)
    )
    if (i == 3L) {
        p <- p + guides(alpha = guide_legend(nrow = 1), size = guide_legend(nrow = 1))
    } else {
        p <- p + guides(colour = guide_colourbar(barwidth = grid::unit(100, 'pt'),
                                                barheight = grid::unit(8, 'pt')),
                        size = guide_legend(nrow = 1))
    }
    assert_data_identity(fig6_block[[i]], p)
    save_figure(paste0('Supplemental_Figure_S10', names6[[i]], '.pdf'), p, 6.15)
}

fig5a <- fig5[[1]][[1]] + readable_theme + labs(tag = NULL) + theme(
    legend.position = 'bottom', legend.box = 'vertical',
    legend.title.position = 'top', legend.box.just = 'left',
    plot.margin = margin(5, 8, 5, 5)
)
assert_data_identity(fig5[[1]][[1]], fig5a)
save_figure('Supplemental_Figure_S10f.pdf', fig5a, 6.6)
writeLines('Data and mappings are identical for all 10 reused component plots; the S10a alpha legend label at rank 35 was corrected from 30.',
           file.path(out, 'plot_data_preservation_audit.txt'))
