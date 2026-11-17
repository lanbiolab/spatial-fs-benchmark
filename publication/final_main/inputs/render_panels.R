suppressPackageStartupMessages({
 library(ggplot2);library(patchwork);library(ggforce);library(colorspace)
})
out <- 'results/revision_annotation_readability_20261109'
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

f3 <- readRDS('results/revision_minor1_3_20261106/figures/figure_num_features_benchmark.rds')
t3 <- theme(axis.text.x=element_text(size=7.8),axis.text.y=element_text(size=7.8),
 axis.title=element_text(size=8.2),strip.text=element_text(size=8.2),
 plot.title=element_text(size=8.6),plot.tag=element_text(size=9.5,face='bold'),
 legend.text=element_text(size=7.8),legend.title=element_text(size=8.1),
 legend.key.width=grid::unit(13,'pt'),legend.key.height=grid::unit(10,'pt'),
 legend.spacing.x=grid::unit(3,'pt'))
p3 <- f3 & t3
p3[[1]] <- p3[[1]] + theme(axis.title.y.left=element_text(vjust=.5,margin=margin(r=5)),
 axis.text.x=element_blank(),strip.text=element_text(size=8.2))
p3[[2]] <- p3[[2]] + theme(axis.text.x=element_blank(),strip.text=element_blank())
p3[[3]] <- p3[[3]] + theme(strip.text=element_blank())
for(i in 1:3) check(f3[[i]],p3[[i]])
# Release only the y-axis label from the common patchwork alignment column.
p3[[1]] <- free(p3[[1]],type='label',side='l')
save_pdf('Figure3.pdf',p3,468,518.4)

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
for(i in 1:2) check(f4[[i]],p4[[i]])
save_pdf('Figure4a.pdf',p4,540,290)

# Original 34-setting biological-context plot, not the later Random-2000 plot.
f5 <- readRDS('results/fig6_biological_context/figure_fig6_biological_context.rds')
p5 <- f5 & theme(axis.text.x=element_text(size=6.2),strip.text=element_text(size=6.7),
 plot.tag=element_text(size=9,face='bold'),plot.margin=margin(2,2,2,2),
 legend.text=element_text(size=6.2),legend.title=element_text(size=6.6),
 legend.key.width=grid::unit(6,'pt'),legend.spacing.x=grid::unit(8,'pt'),
 legend.margin=margin(0,6,0,6),legend.box.spacing=grid::unit(3,'pt'))
for(i in 1:2) {
 p5[[i]][[1]] <- p5[[i]][[1]] + theme(axis.text.y=element_text(size=6.4,colour='black',margin=margin(r=2)))
 for(j in 2:length(f5[[i]])) p5[[i]][[j]] <- p5[[i]][[j]] + theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())
 for(j in 1:length(f5[[i]])) check(f5[[i]][[j]],p5[[i]][[j]])
}
p5[[1]][[2]] <- p5[[1]][[2]] + guides(fill=guide_colourbar(order=2,title.position='top',barwidth=grid::unit(45,'pt'),barheight=grid::unit(6,'pt')))
p5[[1]][[3]] <- p5[[1]][[3]] + guides(colour=guide_colourbar(order=3,title.position='top',barwidth=grid::unit(45,'pt'),barheight=grid::unit(6,'pt')))
for(i in 1:2) p5[[2]][[i]] <- p5[[2]][[i]] + guides(fill=guide_colourbar(order=i,title.position='top',barwidth=grid::unit(55,'pt'),barheight=grid::unit(6,'pt')))
save_pdf('Figure5.pdf',p5,540,612)

f6 <- readRDS('revision_r1_major1_five_20261105/figures/fig6/figure_fig6_spatial_integration.rds')
p6 <- f6
for(i in 1:4) {
 p6[[1]][[i]] <- p6[[1]][[i]] + theme(axis.text.x=element_text(size=7.2),
   strip.text=element_text(size=7.8),plot.tag=element_text(size=9.8,face='bold'),
   plot.margin=margin(2,1,2,1))
 if(i==1) p6[[1]][[i]] <- p6[[1]][[i]] + theme(axis.text.y=element_text(size=7.9,colour='black',margin=margin(r=2)))
 else p6[[1]][[i]] <- p6[[1]][[i]] + theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())
 check(f6[[1]][[i]],p6[[1]][[i]])
}
# Keep the original separate legend grobs and their scale keys unchanged.
save_pdf('Figure6abcd.pdf',p6,590,426)
writeLines('14 ggplot components preserve their data, mappings, layer data and layer mappings. Only themes/layout changed; the existing Figure4 rank-35 legend correction is retained.',file.path(out,'data_preservation.txt'))
