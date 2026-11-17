suppressPackageStartupMessages(library(ggplot2))
p <- readRDS("revision_r1_major1_20261105/figures_generated/original_script_fig3/figure_num_features_benchmark.rds")
out <- "revision_r1_major1_20261105/figures_generated/original_script_fig3/figure_num_features_benchmark"
ggsave(paste0(out,".svg"),p,width=7.8,height=6.8,device=svglite::svglite,bg="white")
ggsave(paste0(out,".tiff"),p,width=7.8,height=6.8,dpi=600,device=ragg::agg_tiff,compression="lzw",bg="white")
