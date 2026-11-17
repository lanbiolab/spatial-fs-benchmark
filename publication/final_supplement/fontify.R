suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(ggforce);library(colorspace);library(grid);library(cowplot)})
out='results/revision_all_typography_20261109'
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
