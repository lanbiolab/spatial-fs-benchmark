#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(ggplot2);library(dplyr);library(readr);library(ggridges);library(patchwork)})
args <- commandArgs(trailingOnly=TRUE);dir <- args[[1]]
lev <- c('bASW','dASW','iLISI','dLISI','ILL','GC','ari','nmi','CHAOS','PAS','silhouette')
labels <- c(ari='ARI',nmi='NMI',silhouette='Silhouette')
v <- read_tsv(file.path(dir,'metric_characterization_values.tsv'),show_col_types=FALSE) |>
 mutate(metric=factor(.data$metric_name,levels=rev(lev))) |>
 group_by(.data$metric) |> mutate(metric_mean=mean(.data$ScaledValue)) |> ungroup()
f <- read_tsv(file.path(dir,'feature_budget_correlations.tsv'),show_col_types=FALSE) |>
 mutate(metric=factor(.data$metric,levels=rev(lev)))
summ <- f |> group_by(.data$metric) |> summarise(mean_rho=mean(.data$rho),.groups='drop')
c <- read_tsv(file.path(dir,'metric_pair_correlations_summary.tsv'),show_col_types=FALSE) |>
 mutate(metric_a=factor(.data$metric_a,levels=rev(lev)),metric_b=factor(.data$metric_b,levels=lev))
stopifnot(nrow(c)==121,all(c$n==6),nrow(f)==66)
base <- theme_classic(base_size=8,base_family='Helvetica')+theme(plot.title=element_text(size=8,face='bold',hjust=.5),axis.title=element_text(size=7),axis.text=element_text(size=6.5,colour='black'),plot.margin=margin(2,2,2,2))
display <- function(x) ifelse(x %in% names(labels),unname(labels[x]),x)
p1 <- ggplot(v,aes(x=.data$ScaledValue,y=.data$metric,group=.data$metric,fill=.data$metric_mean))+
 geom_density_ridges(scale=.86,rel_min_height=.01,colour='white',linewidth=.2,bandwidth=.05,from=0,to=1)+
 scale_fill_viridis_c(option='magma',limits=c(0,1),guide='none')+
 scale_y_discrete(labels=display)+scale_x_continuous(limits=c(0,1),breaks=c(0,.5,1))+
 labs(title='Observed\nscore distribution',x='Frozen scaled value',y=NULL)+base+
 theme(axis.text.y=element_text(size=6.6),axis.ticks.y=element_blank())
p2 <- ggplot(f,aes(x=.data$rho,y=.data$metric))+
 geom_vline(xintercept=0,colour='grey70',linewidth=.3)+
 geom_point(position=position_jitter(width=0,height=.09,seed=27),size=.9,colour='#659ebc',alpha=.8)+
 geom_point(data=summ,aes(x=.data$mean_rho,y=.data$metric),shape=23,fill='#e69f00',colour='white',size=1.8,stroke=.25)+
 scale_x_continuous(limits=c(-1,1),breaks=c(-1,0,1))+scale_y_discrete(labels=display)+
 labs(title='Correlation with\nfeature budget',x='Spearman correlation',y=NULL)+base+
 theme(axis.text.y=element_blank(),axis.ticks.y=element_blank())
p3 <- ggplot(c,aes(x=.data$metric_b,y=.data$metric_a,fill=.data$mean_rho))+
 geom_tile(colour='white',linewidth=.2)+scale_fill_gradient2(low='#398363',mid='white',high='#77408d',midpoint=0,limits=c(-1,1),name='Mean\ncorrelation')+
 scale_x_discrete(labels=display)+scale_y_discrete(labels=display)+coord_fixed()+
 labs(title='Correlation between metrics',x=NULL,y=NULL)+base+
 theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.text.x=element_text(angle=90,hjust=1,vjust=.5,size=6),legend.position='bottom',legend.title=element_text(size=6),legend.text=element_text(size=6),legend.key.width=unit(.6,'cm'))
fig <- (p1|p2|p3)+plot_layout(widths=c(1.35,1.15,3.3),guides='collect') & theme(legend.position='bottom')
saveRDS(fig,file.path(dir,'metric_corrected.rds'))
ggsave(file.path(dir,'metric_corrected.png'),fig,width=11.2,height=6.2,dpi=300,bg='white')
ggsave(file.path(dir,'metric_corrected.pdf'),fig,width=11.2,height=6.2,device=cairo_pdf,bg='white')
