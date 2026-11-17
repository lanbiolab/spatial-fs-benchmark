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
save(p1b,'Figure1b',360)

f2=readRDS('results/manuscript_rebuild_v2/main_figure2_validation/Figure2_spatial_validation.rds')
p2=resize(f2,6.5,6.5)
p2[[1]][[1]]=p2[[1]][[1]]+theme(panel.spacing.x=unit(12,'pt'))
p2[[1]][[2]]=p2[[1]][[2]]+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),strip.text=element_text(size=5.5))
p2[[2]][[2]]$facet$params$labeller=as_labeller(c(Jaccard='Top-N\nJaccard',MeanReferencePercentile='Percentile',RankSpearman='Rank\nagreement'))
p2[[2]][[2]]=p2[[2]][[2]]+labs(x="Rank using Moran's I\nreference")
p2[[2]][[2]]=p2[[2]][[2]]+theme(axis.text.x=element_text(size=5.8),strip.text=element_text(size=6),panel.spacing.x=unit(4,'pt'))
p2[[2]][[3]]=p2[[2]][[3]]+labs(x='Association of\nreproducibility\nwith ARI')
for(j in seq_along(p2[[2]][[2]]$layers))if(inherits(p2[[2]][[2]]$layers[[j]]$geom,'GeomText'))p2[[2]][[2]]$layers[[j]]$aes_params$size=5.2/ggplot2::.pt
save(p2,'Figure2',565)

# Remaining blocks of Figure 4 retain their source data, scales and plot types.
f4b=readRDS('results/reference_aligned_v2/fig4bcd/figure_fig4b_panel_b.rds')
f4c=readRDS('results/reference_aligned_v2/fig4bcd/figure_fig4c_panel_c.rds')
f4d=readRDS('results/reference_aligned_v2/fig4bcd/figure_fig4d_panel_d.rds')
for(n in c('b','c','d')) {
 p=resize(get(paste0('f4',n)),5.6,5.5)+labs(tag=n,x=NULL,y=NULL)+theme(plot.tag=element_text(size=9,face='bold'),plot.tag.position=c(.005,.995))
 if(n=='b')p=p+theme(legend.box='horizontal',legend.key.width=unit(5,'pt'),legend.spacing.x=unit(3,'pt'),legend.margin=margin(0,1,0,1))
 if(n=='c')p=p+labs(x='Number of selected\nfeatures (log scale)')+theme(legend.position='bottom',legend.direction='horizontal',legend.title.position='top',legend.key.width=unit(5,'pt'),legend.margin=margin(0,0,0,0))+guides(colour=guide_legend(nrow=1))
 save(p,paste0('Figure4',n),210,c(b=252,c=126,d=162)[n])
}
f4e=readRDS('results/revision_major2_response_20261105/figures/fig4root/fig4e_hybrid/figure_fig4e_hybrid_panel.rds')
p4e=resize(f4e,5.8,6)
p4e[[1]]=p4e[[1]]+labs(tag='e')+theme(plot.tag=element_text(size=9,face='bold'))
save(p4e,'Figure4e',136)

# Figure 6 legends are rebuilt from the unchanged source scales and aligned
# below their corresponding data blocks, excluding the method-label column.
f6=readRDS('revision_r1_major1_five_20261105/figures/fig6/figure_fig6_spatial_integration.rds')
top=f6[[1]]
for(i in 1:4){top[[i]]=resize(top[[i]],7.2)+theme(strip.text=element_text(size=7.8),plot.tag=element_text(size=9.8,face='bold'),plot.margin=margin(2,1,2,1));if(i==1)top[[i]]=top[[i]]+theme(axis.text.y=element_text(size=7.9,margin=margin(r=2)))}
lt=theme(legend.position='bottom',legend.title.position='top',legend.box='horizontal',legend.text=element_text(size=6.8),legend.title=element_text(size=7.2,face='bold'),legend.margin=margin(0,1,0,1),legend.box.margin=margin(0,0,0,0),legend.key.width=unit(8,'pt'),legend.key.height=unit(7,'pt'))
upper=lower=list()
titles=c('Mean score','Mean difference\nto scVI','Mean rank','Rank difference\nto scVI')
for(i in 1:4){
 z=f6[[1]][[i]]+lt
 g=if(i==3)guides(alpha=guide_legend(title=titles[i],nrow=2,byrow=TRUE),size='none') else guides(colour=guide_colourbar(title=titles[i],barwidth=unit(if(i %in% c(2,4))55 else 65,'pt'),barheight=unit(6,'pt')),size='none')
 upper[[i]]=wrap_elements(cowplot::get_legend(z+g))
 lower[[i]]=wrap_elements(cowplot::get_legend(z+theme(legend.key.width=unit(5,'pt'),legend.text=element_text(size=6.5))+guides(colour='none',alpha='none',size=guide_legend(title=if(i==3)'Rank s.d.' else 's.d.',nrow=1))))
}
row1=wrap_plots(c(list(plot_spacer()),upper,list(plot_spacer())),nrow=1,widths=c(149,125,90,136,80,10))
row2=wrap_plots(c(list(plot_spacer()),lower,list(plot_spacer())),nrow=1,widths=c(149,125,90,136,80,10))
p6=wrap_plots(top,row1,row2,ncol=1,heights=c(1,.085,.075))
save(p6,'Figure6abcd',438,590)

# S1 is the archived seed-stability plot (historical build name S3).
save(resize(readRDS('results/manuscript_rebuild_v2/supplemental/figures/Supplemental_Fig_S3.rds'),7),'Supplemental_Figure_S1',510)

# Reconstruct S2 from the already-computed five-study rank tables, stopping before exports.
env=new.env(parent=globalenv());env$commandArgs=function(...)c('revision_r1_major1_five_20261105/analysis/robustness_inputs',out)
code=readLines('revision_r1_major1_five_20261105/analysis/plot_s2_five.R')
eval(parse(text=code[seq_len(grep('^width_mm <-',code)[1]-1)]),env)
save(resize(env$figure,7),'Supplemental_Figure_S2',560)

# S3 is split into a and b pages: no panel or legend is removed.
f3=readRDS('results/reference_aligned_v2/baselines/figure_baselines.rds')
p3a=resize(attr(f3[[1]],'grobs')$full,7)+labs(tag='a')+theme(plot.tag=element_text(size=9,face='bold'),legend.position='bottom')
p3a$facet$params$ncol=3;p3a$facet$params$nrow=NULL
save(p3a,'Supplemental_Figure_S3a',565)
p3b=resize(attr(f3[[2]],'grobs')$full,6.8)
p3b[[1]]=p3b[[1]]+labs(tag='b')+theme(panel.spacing.x=unit(5,'pt'),axis.text.x=element_text(size=5.8))
p3b[[2]]=p3b[[2]]+labs(title='2  Scale using frozen\nobserved ranges')
p3b[[3]]=p3b[[3]]+labs(title='3  Average by\nmetric type')
p3b[[4]]=p3b[[4]]+labs(title='4  CoreOverall\nscore')
p3b=p3b+plot_annotation(caption=paste(strwrap(f3$patches$annotation$caption,width=130),collapse='\n')) & theme(plot.caption=element_text(size=6.5,hjust=0))
save(p3b,'Supplemental_Figure_S3b',420)

save(resize(readRDS('revision_r1_major1_five_20261105/figures/batch_controls_seven/batch_mixing_controls_seven.rds'),7.5),'Supplemental_Figure_S5',440)
save(resize(readRDS('revision_r1_major1_five_20261105/figures/embryo_label_free/embryo_label_free_sensitivity.rds'),7.3),'Supplemental_Figure_S6',560)
for(n in c(a='heldout_ARI',b='external_markers')) {
 p=readRDS(paste0('results/revision_major4_response_20261106/figures/Supplemental_Major4_',n,'.rds'))
 save(resize(p,7.4),paste0('Supplemental_Figure_S7',names(which(c(a='heldout_ARI',b='external_markers')==n))),355)
}

# Existing plotting scripts read saved result tables; export destinations are intercepted.
capture_script=function(path,handler) {
 e=new.env(parent=globalenv());e$ggsave=function(filename,plot,...){if(grepl('\\.pdf$',filename))handler(plot,basename(filename))}
 # Redirect any source-data export into this isolated run too.
 e$write_tsv=function(x,file,...)invisible(NULL)
 e$writeLines=function(text,con,...)invisible(NULL)
 eval(parse(text=readLines(path)),e)
}
capture_script('results/revision_r3_comment5_20261107/scripts/plot_results.R',function(p,name){
 mapping=c(R3C5_common_ranks_scvi='S9a',R3C5_common_ranks_cellcharter='S9b',R3C5_common_ranks_graphst='S9c',R3C5_annotation_ARI='S9d',R3C5_platform_correlations='S9e')
 key=mapping[sub('\\.pdf$','',name)];fs=if(key=='S9d')6.7 else 7.4
 save(resize(p,fs),paste0('Supplemental_Figure_',key),if(key=='S9d')612 else if(key=='S9e')350 else 500)
})
capture_script(file.path(out,'plot_platforms_original.R'),function(p,name){
 key=if(grepl('common_component',name))'S4a' else 'S4b'
 save(resize(p,7.4),paste0('Supplemental_Figure_',key),if(key=='S4a')430 else 300)
})
capture_script(file.path(out,'plot_individual_original.R'),function(p,name){
 key=if(grepl('_ARI',name))'S8a' else 'S8b'
 save(resize(p,7.4),paste0('Supplemental_Figure_',key),if(key=='S8a')370 else 540)
})
capture_script(file.path(out,'plot_targets_original.R'),function(p,name){save(resize(p,7.4),'Supplemental_Figure_S12',530)})
code=readLines('results/revision_minor1_3_20261106/scripts/plot_minor3_original_style.R')
start=grep('^raw <- readr::read_tsv',code)[1];end=grep('^save_figure_files\\(budget_plot',code)[1]-1
eval(parse(text=code[start:end]))
save(resize(budget_plot,7.4),'Supplemental_Figure_S11',300)
writeLines(audit,file.path(out,'data_preservation.txt'))
