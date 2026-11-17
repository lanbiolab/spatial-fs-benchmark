suppressPackageStartupMessages({
 library(ggplot2);library(patchwork);library(ggforce);library(colorspace)
})
out <- 'results/revision_figure4_rank_spacing_20261109'
dir.create(out, recursive=TRUE, showWarnings=FALSE)
save_pdf <- function(name,p,w,h) {
 set.seed(20261109)
 ggsave(file.path(out,name),p,width=w/72,height=h/72,units='in',
        device=grDevices::cairo_pdf,family='Arial',bg='white')
}
check <- function(a,b) {
 stopifnot(identical(a$data,b$data),identical(a$mapping,b$mapping),length(a$layers)==length(b$layers))
 for(j in seq_along(a$layers)) stopifnot(
   identical(a$layers[[j]]$layer_data(a$data),b$layers[[j]]$layer_data(b$data)),
   identical(a$layers[[j]]$mapping,b$layers[[j]]$mapping))
}

f4 <- readRDS('results/revision_major2_response_20261105/figures/fig4a/figure_fig4a_spatial_benchmark.rds')
p4 <- f4 & theme(axis.text.x=element_text(size=6.2),strip.text=element_text(size=7),
 plot.title=element_text(size=7.2),legend.text=element_text(size=6.2),
 legend.title=element_text(size=6.6),legend.key.width=grid::unit(10,'pt'),
 legend.key.height=grid::unit(8,'pt'),legend.spacing.x=grid::unit(3,'pt'),
 legend.box.spacing=grid::unit(2,'pt'),plot.margin=margin(1,1,1,1))
p4[[1]] <- p4[[1]] + labs(tag='a') + theme(axis.text.y=element_text(size=6.4,colour='black',margin=margin(r=2)),
 plot.tag=element_text(size=9,face='bold',vjust=1),plot.tag.position=c(.006,.994))
p4[[2]] <- p4[[2]] + theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())
alpha <- f4[[2]]$scales$get_scales('alpha')$clone()
alpha$labels <- as.character(alpha$breaks)
p4[[2]] <- p4[[2]] + alpha
size_scale <- f4[[2]]$scales$get_scales('size')$clone()
original_palette <- size_scale$palette
size_scale$palette <- function(x) original_palette(x)*0.55
p4[[2]] <- p4[[2]] + size_scale
for(i in 1:2) check(f4[[i]],p4[[i]])
stopifnot(all.equal(p4[[2]]$scales$get_scales('size')$palette(c(0,.25,.5,.75,1)),f4[[2]]$scales$get_scales('size')$palette(c(0,.25,.5,.75,1))*.55))
saveRDS(p4,file.path(out,'Figure4a.rds'))
writeLines('Data, mapping, layer data and ordering unchanged. SDRank size palette multiplied by 0.55, including its automatic legend. MeanRank alpha and Type colour unchanged.',file.path(out,'verification.txt'))
save_pdf('Figure4a.pdf',p4,540,290)

