original <- readLines('scripts/generate_main_spatial_validation_figure.R')
eval(parse(text=paste(original[seq_len(grep('^semi <-',original)[1]-1L)],collapse='\n')))
library(patchwork)
root <- 'results/revision_r2_m1_response_20261111'
a <- file.path(root,'analysis'); out <- file.path(root,'figures'); dir.create(out,recursive=TRUE,showWarnings=FALSE)
theme_scale <- base_theme + theme(text=element_text(family='Arial'),axis.text=element_text(size=7.5),axis.title=element_text(size=8),plot.title=element_text(size=8.5,face='bold'),plot.tag=element_text(size=10,face='bold'),plot.margin=margin(2,2,2,2,'mm'))
fine <- read_tsv(file.path(a,'geometry_fine_footprint.tsv'),show_col_types=FALSE)
coarse <- read_tsv(file.path(a,'geometry_pseudo_visium.tsv'),show_col_types=FALSE)
circles <- bind_rows(lapply(seq_len(nrow(coarse)),function(i){t<-seq(0,2*pi,length.out=101);data.frame(x=coarse$x[i]+27.5*cos(t),y=coarse$y[i]+27.5*sin(t),id=i)}))
geo <- function(p) p+coord_equal(xlim=c(-220,220),ylim=c(-220,220),expand=FALSE)+scale_x_continuous(breaks=c(-200,0,200))+scale_y_continuous(breaks=c(-200,0,200))+labs(x='Relative x (nominal \u00b5m)',y='Relative y (nominal \u00b5m)')+theme_scale
p <- geo(ggplot()+geom_path(data=circles,aes(x,y,group=id),colour='grey65',linewidth=.25)+geom_point(data=fine,aes(x,y),size=.65,colour='#e78ac3'))+labs(title='Same beads, unpooled')
q <- geo(ggplot()+geom_polygon(data=circles,aes(x,y,group=id),fill='#377eb8',colour=NA,alpha=.85))+labs(title='55-\u00b5m pseudo-spots')
u <- read_tsv(file.path(a,'aggregation_audit.tsv'),show_col_types=FALSE) |> select(sample_id,fine_footprint_median_knn8_um_assumed,pseudo_visium_median_knn8_um_assumed) |> pivot_longer(-sample_id,names_to='arm',values_to='distance') |> mutate(arm=factor(arm,levels=c('fine_footprint_median_knn8_um_assumed','pseudo_visium_median_knn8_um_assumed'),labels=c('Fine','Pooled')))
r <- ggplot(u,aes(arm,distance,group=sample_id))+geom_hline(yintercept=100*sqrt(3),linetype='dashed',colour='grey45',linewidth=.35)+geom_line(colour='grey65',linewidth=.35)+geom_point(aes(colour=arm),size=1.8,position=position_jitter(width=.045,height=0,seed=11))+scale_colour_manual(values=c(Fine='#e78ac3',Pooled='#377eb8'),guide='none')+scale_y_continuous(limits=c(0,200),breaks=c(0,50,100,150,200))+labs(x=NULL,y='Median 8th-neighbour distance\n(nominal \u00b5m)',title='Neighbourhood scale')+theme_scale
d <- read_tsv(file.path(a,'paired_rank_correlations.tsv'),show_col_types=FALSE)
lev <- c('paired_complete','exact_2000_both','paired_excluding_nnsvg')
met <- c('PAS','CHAOS','silhouette','iLISI','bASW')
d <- d |> mutate(panel=factor(panel,levels=lev),workflow=factor(workflow,levels=c('scvi','cellcharter'),labels=c('scVI','CellCharter')),y=match(metric,met)+c(.20,0,-.20)[match(panel,lev)])
s <- ggplot(d,aes(rank_rho,y,colour=panel,shape=panel))+geom_vline(xintercept=0,colour='grey55',linetype='dashed',linewidth=.3)+geom_point(size=2)+geom_text(data=d |> filter(panel=='paired_complete'),aes(label=sprintf('%.2f',rank_rho)),nudge_x=.08,hjust=0,size=2.4,colour='grey15',show.legend=FALSE)+facet_wrap(~workflow,nrow=1)+scale_x_continuous(limits=c(-1,1),breaks=c(-1,-.5,0,.5,1))+scale_y_continuous(breaks=seq_along(met),labels=c('PAS','CHAOS','Silhouette','iLISI','bASW'),limits=c(.55,5.45))+scale_colour_manual(values=c('#377eb8','#e78ac3','#555555'),labels=c('Paired methods (n = 28)','Exact 2,000 (n = 26)','Excluding nnSVG (n = 27)'),name=NULL)+scale_shape_manual(values=c(16,17,15),labels=c('Paired methods (n = 28)','Exact 2,000 (n = 26)','Excluding nnSVG (n = 27)'),name=NULL)+labs(x='Fine-pooled rank correlation (Spearman rho)',y=NULL)+theme_scale+theme(legend.position='bottom',legend.text=element_text(size=7.5),legend.key.width=unit(3,'mm'),legend.spacing.x=unit(1,'mm'),panel.grid.major.y=element_line(colour='grey93',linewidth=.25),panel.grid.major.x=element_blank(),panel.grid.minor=element_blank(),strip.text=element_text(size=8.5,face='bold'),panel.spacing=unit(7,'mm'))
fig <- ((p+q+r)+plot_layout(widths=c(1,1,.9)))/s+plot_layout(heights=c(1,.82))+plot_annotation(tag_levels='a')
for(ext in c('pdf','png','svg')) {
 device <- switch(ext,pdf=grDevices::cairo_pdf,png=ragg::agg_png,svg=svglite::svglite)
 args <- list(filename=file.path(out,paste0('R2_M1_scale_design_summary.',ext)),plot=fig,width=183,height=145,units='mm',device=device,bg='white')
 if(ext=='pdf')args$family<-'Arial'
 if(ext=='png')args$dpi<-300
 do.call(ggsave,args)
}
