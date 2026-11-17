suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(ggforce);library(colorspace);library(dplyr);library(readr);library(cowplot);library(grid);library(ggridges)})
out <- 'results/revision_all_figure_readability_20261109'
dir.create(out,recursive=TRUE,showWarnings=FALSE)
audit <- character()
save <- function(p,name,h,w=540) {
 message('Saving ',name)
 set.seed(20261109)
 ggsave(file.path(out,paste0(name,'.pdf')),p,width=w/72,height=h/72,units='in',device=grDevices::cairo_pdf,family='Arial',bg='white')
 saveRDS(p,file.path(out,paste0(name,'.rds')))
}
check <- function(a,b) {
 expr=function(x)lapply(x,function(y)paste(deparse(rlang::get_expr(y)),collapse=''))
 stopifnot(identical(a$data,b$data),identical(expr(a$mapping),expr(b$mapping)),length(a$layers)==length(b$layers))
 for(j in seq_along(a$layers))stopifnot(identical(a$layers[[j]]$layer_data(a$data),b$layers[[j]]$layer_data(b$data)),identical(expr(a$layers[[j]]$mapping),expr(b$layers[[j]]$mapping)))
 audit <<- c(audit,paste('Preserved data/mappings/layers:',nrow(a$data),length(a$layers)))
}
grob_fonts <- function(g,fs) {
 if(inherits(g,'text')) {
  if(identical(g$label,'Spatial transcriptomics\ndatasets'))g$label='Spatial\ntranscriptomics\ndatasets'
  if(identical(g$label,'Expression-driven / HVG'))g$label='Expression-driven\n/ HVG'
  if(identical(g$label,'Spatially informed / SVG'))g$label='Spatially informed\n/ SVG'
  if(identical(g$label,'Cross-statistic'))g$label='Cross-\nstatistic'
  if(identical(g$label,'Feature-number\neffects'))g$label='Feature-\nnumber\neffects'
  if(identical(g$label,'Integration-method\ncomparison'))g$label='Integration-\nmethod\ncomparison'
  g$label=gsub('Mean feature-number correlation','Mean feature-number\ncorrelation',g$label,fixed=TRUE)
  g$label=gsub('Mean dataset-feature correlation','Mean dataset-feature\ncorrelation',g$label,fixed=TRUE)
 }
 if(inherits(g,'text') && !is.null(g$gp$fontsize))g$gp$fontsize=pmax(g$gp$fontsize,fs)
 if(!is.null(g$grobs))g$grobs=lapply(g$grobs,grob_fonts,fs=fs)
 if(!is.null(g$children))for(i in seq_along(g$children))g$children[[i]]=grob_fonts(g$children[[i]],fs)
 g
}
resize <- function(p,fs=7,legend=fs) {
 old=p
 if(inherits(p,'patchwork')) {for(i in seq_len(length(p)))p[[i]]=resize(p[[i]],fs,legend);return(p)}
 if(inherits(p,'wrapped_patch')) {
  gs=attr(p,'grobs');for(n in names(gs))if(inherits(gs[[n]],'ggplot'))gs[[n]]=resize(gs[[n]],fs,legend) else if(inherits(gs[[n]],'grob'))gs[[n]]=grob_fonts(gs[[n]],fs)
  attr(p,'grobs')=gs;return(p)
 }
 if(!inherits(p,'ggplot'))return(p)
 p=unserialize(serialize(p,NULL))
 for(k in c('axis.text.x','axis.text.y','axis.title.x','axis.title.y','strip.text','plot.subtitle','plot.title','legend.text','legend.title')) {
  if(!inherits(ggplot2::calc_element(k,p$theme),'element_blank'))p=p+do.call(theme,setNames(list(element_text(size=if(grepl('legend',k))legend else if(grepl('title',k))fs+.5 else fs)),k))
 }
 p=p+theme(plot.margin=margin(3,3,3,3),legend.margin=margin(2,4,2,4),legend.spacing.x=unit(7,'pt'),legend.box.spacing=unit(4,'pt'))
 for(j in seq_along(p$layers)) {
  if(inherits(p$layers[[j]]$geom,'GeomText') && !is.null(p$layers[[j]]$aes_params$size))p$layers[[j]]$aes_params$size=pmax(p$layers[[j]]$aes_params$size,fs/ggplot2::.pt)
  if(inherits(p$layers[[j]]$geom,'GeomText') && !is.null(p$layers[[j]]$aes_params$label)) {
   p$layers[[j]]$aes_params$label=gsub('Mean feature-number correlation','Mean feature-number\ncorrelation',p$layers[[j]]$aes_params$label,fixed=TRUE)
   p$layers[[j]]$aes_params$label=gsub('Mean dataset-feature correlation','Mean dataset-feature\ncorrelation',p$layers[[j]]$aes_params$label,fixed=TRUE)
  }
  if(!is.null(p$layers[[j]]$geom_params$grob))p$layers[[j]]$geom_params$grob=grob_fonts(p$layers[[j]]$geom_params$grob,fs)
 }
 check(old,p);p
}

# Main 1: exact original diagram grob and original metric plot objects.
f1=readRDS('results/figure1_complete_v6/Figure1_panels.rds')
save(grob_fonts(f1$dag,6.0),'Figure1a',220)
f1b=readRDS('results/reference_aligned_v2/metric_overview/figure_metric_overview.rds')
p1b=resize(f1b,6.4,6.2)
for(i in seq_len(length(p1b[[2]])))for(j in seq_along(p1b[[2]][[i]]$layers)) {
 l=p1b[[2]][[i]]$layers[[j]]
 if(inherits(l$geom,'GeomText') && length(l$aes_params$label)>1 && all(grepl('^-?[0-9.]+$',l$aes_params$label)))p1b[[2]][[i]]$layers[[j]]$aes_params$size=5.8/ggplot2::.pt
}
for(i in 1:3)p1b[[1]][[i]]=p1b[[1]][[i]]+theme(plot.margin=margin(3,6,3,5))
p1b[[1]][[4]]=p1b[[1]][[4]]+theme(strip.text.x=element_text(size=5.5),strip.text.y=element_text(size=6.4))
save(p1b,'Figure1b',360)

