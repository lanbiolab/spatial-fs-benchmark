suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(ggforce);library(colorspace);library(grid);library(cowplot)})
out='results/revision_arial_spacing_20261109'
options(warn=1)
audit=character()
fontify_grob=function(g){
 if(inherits(g,'text'))g$gp$fontfamily='Arial'
 if(!is.null(g$grobs))g$grobs=lapply(g$grobs,fontify_grob)
 if(!is.null(g$children))for(i in seq_along(g$children))g$children[[i]]=fontify_grob(g$children[[i]])
 g
}
fontify=function(p){
 if(inherits(p,'grob'))return(fontify_grob(p))
 if(inherits(p,'wrapped_patch')){
  gs=attr(p,'grobs');for(n in names(gs))gs[[n]]=fontify(gs[[n]])
  attr(p,'grobs')=gs
 }
 if(inherits(p,'patchwork')){
  for(i in seq_len(length(p)))p[[i]]=fontify(p[[i]])
 }
 if(!inherits(p,'ggplot'))return(p)
 old=unserialize(serialize(p,NULL));p=unserialize(serialize(p,NULL))
 for(n in names(p$theme))if(inherits(p$theme[[n]],'element_text'))p$theme[[n]]$family='Arial'
 p=p+theme(text=element_text(family='Arial'))
 for(i in seq_along(p$layers)){
  if(inherits(p$layers[[i]]$geom,'GeomText')||inherits(p$layers[[i]]$geom,'GeomLabel'))p$layers[[i]]$aes_params$family='Arial'
  if(!is.null(p$layers[[i]]$geom_params$grob))p$layers[[i]]$geom_params$grob=fontify_grob(p$layers[[i]]$geom_params$grob)
  stopifnot(identical(old$layers[[i]]$layer_data(old$data),p$layers[[i]]$layer_data(p$data)))
 }
 mapping_expr=function(x)lapply(x,function(y)paste(deparse(rlang::get_expr(y)),collapse=''))
 stopifnot(identical(old$data,p$data),identical(mapping_expr(old$mapping),mapping_expr(p$mapping)))
 audit<<-c(audit,'Data and mappings unchanged; font family normalized to Arial.')
 p
}
emit=function(p,name,w,h){
 set.seed(20261109)
 ggplot2::ggsave(file.path(out,name),fontify(p),width=w/72,height=h/72,units='in',device=grDevices::cairo_pdf,family='Arial',bg='white')
}
emit(readRDS('results/revision_all_figure_readability_20261109/Figure1a.rds'),'Figure1a.pdf',540,220)
emit(readRDS('results/revision_all_figure_readability_20261109/Figure1b.rds'),'Figure1b.pdf',540,360)
original_palette=readRDS('results/revision_major2_response_20261105/figures/fig4a/figure_fig4a_spatial_benchmark.rds')[[2]]$scales$get_scales('size')$palette
emit(readRDS('results/revision_figure4_rank_spacing_20261109/Figure4a.rds'),'Figure4a.pdf',540,290)
emit(readRDS('results/revision_all_figure_readability_20261109/Figure6abcd.rds'),'Figure6abcd.pdf',590,438)
emit(readRDS('results/revision_r3_comment6_20261108/figures/Figure6e_task_guide.rds'),'Figure6e.pdf',590,154)
capture=function(path,keep){
 e=new.env(parent=globalenv())
 e$ggsave=function(filename,plot,...){
  n=basename(filename);if(!n %in% keep)return(invisible(NULL))
  args=list(...);args$filename=file.path(out,n);args$plot=fontify(plot)
  args$device=grDevices::cairo_pdf;args$family='Arial';args$bg='white'
  do.call(ggplot2::ggsave,args)
 }
 e$writeLines=function(...)invisible(NULL)
 eval(parse(text=readLines(path)),e)
}
capture(file.path(out,'inputs/render_panels.R'),c('Figure3.pdf','Figure5.pdf'))
capture(file.path(out,'inputs/render_four_items.R'),'Figure4d.pdf')
capture(file.path(out,'inputs/render_recheck.R'),c('Figure2.pdf','Figure4b.pdf','Figure4c.pdf','Figure4e.pdf'))
writeLines(audit,file.path(out,'font_data_verification.txt'))
