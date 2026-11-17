original <- readLines('scripts/generate_main_spatial_validation_figure.R')
eval(parse(text=paste(original[seq_len(grep('^semi <-',original)[1]-1L)],collapse='\n')))
library(patchwork)
root <- 'results/revision_r2_m1_scale_match_20261111'
a <- file.path(root,'analysis'); out <- file.path(root,'figures'); dir.create(out,recursive=TRUE,showWarnings=FALSE)
save_panel <- function(p,name,h) {
  ggsave(file.path(out,paste0(name,'.pdf')),p,width=183,height=h,units='mm',device=grDevices::cairo_pdf,family='Arial',bg='white')
  ggsave(file.path(out,paste0(name,'.png')),p,width=183,height=h,units='mm',dpi=300,device=ragg::agg_png,bg='white')
  ggsave(file.path(out,paste0(name,'.svg')),p,width=183,height=h,units='mm',device=svglite::svglite,bg='white')
}
theme_scale <- base_theme + theme(text=element_text(family='Arial'),axis.text=element_text(size=7.5),axis.title=element_text(size=8),plot.title=element_text(size=8.5,face='bold'),plot.tag=element_text(size=10,face='bold'),plot.margin=margin(2,2,2,2,'mm'))
fine <- read_tsv(file.path(a,'geometry_fine_footprint.tsv'),show_col_types=FALSE)
coarse <- read_tsv(file.path(a,'geometry_pseudo_visium.tsv'),show_col_types=FALSE)
circles <- bind_rows(lapply(seq_len(nrow(coarse)),function(i){t<-seq(0,2*pi,length.out=101);data.frame(x=coarse$x[i]+27.5*cos(t),y=coarse$y[i]+27.5*sin(t),id=i)}))
geometry <- function(p) p+coord_equal(xlim=c(-220,220),ylim=c(-220,220),expand=FALSE)+scale_x_continuous(breaks=c(-200,0,200))+scale_y_continuous(breaks=c(-200,0,200))+labs(x='Relative x (assumed \u00b5m)',y='Relative y (assumed \u00b5m)')+theme_scale
p <- geometry(ggplot()+geom_path(data=circles,aes(x,y,group=id),colour='grey65',linewidth=.25)+geom_point(data=fine,aes(x,y),size=.65,colour='#e78ac3'))+labs(title='Same beads, unpooled')
q <- geometry(ggplot()+geom_polygon(data=circles,aes(x,y,group=id),fill='#377eb8',colour=NA,alpha=.85))+labs(title='55-\u00b5m pseudo-spots')
u <- read_tsv(file.path(a,'aggregation_audit.tsv'),show_col_types=FALSE) |> select(sample_id,fine_footprint_median_knn8_um_assumed,pseudo_visium_median_knn8_um_assumed) |> pivot_longer(-sample_id,names_to='arm',values_to='distance') |> mutate(arm=factor(arm,levels=c('fine_footprint_median_knn8_um_assumed','pseudo_visium_median_knn8_um_assumed'),labels=c('Fine','Pooled')))
r <- ggplot(u,aes(arm,distance,group=sample_id))+geom_hline(yintercept=100*sqrt(3),linetype='dashed',colour='grey45',linewidth=.35)+geom_line(colour='grey65',linewidth=.35)+geom_point(aes(colour=arm),size=1.8,position=position_jitter(width=.045,height=0,seed=11))+scale_colour_manual(values=c(Fine='#e78ac3',Pooled='#377eb8'),guide='none')+scale_y_continuous(limits=c(0,200),breaks=c(0,50,100,150,200))+labs(x=NULL,y='Median 8th-neighbour distance\n(assumed \u00b5m)',title='Neighbourhood scale')+theme_scale
save_panel((p+q+r)+plot_layout(widths=c(1,1,.9))+plot_annotation(tag_levels='a'),'Scale_Matched_Geometry',83)
if(file.exists(file.path(a,'endpoint_ranks.tsv'))) {
 d <- read_tsv(file.path(a,'endpoint_ranks.tsv'),show_col_types=FALSE)
 nr <- n_distinct(d$fs_method)
 order <- sort(unique(d$fs_method));labs <- method_label(order)
 d <- d |> mutate(Method=factor(fs_method,levels=rev(order),labels=rev(labs)),Arm=factor(arm,levels=c('fine_footprint','pseudo_visium'),labels=c('Fine','Pooled')),Metric=factor(metric_name,levels=c('bASW','iLISI','silhouette','CHAOS','PAS'),labels=c('bASW','iLISI','Silhouette','CHAOS','PAS')))
 for(flow in c('scvi','cellcharter')) {
  x<-d |> filter(integration_method==flow)
  p<-ggplot(x,aes(Arm,Method,fill=rank))+geom_tile(colour=NA)+geom_point(data=x |> filter(rank<=3),shape=21,size=1.1,stroke=.2,fill='white',colour='grey20')+facet_grid(.~Metric)+scale_fill_gradientn(colours=c('#2166ac','#67a9cf','#d1e5f0','#f7f7f7'),limits=c(1,nr),breaks=unique(c(1,10,20,nr)),name='Rank (1 = best)')+labs(x=NULL,y=NULL)+theme_scale+theme(axis.text.y=element_text(size=7.5),axis.text.x=element_text(size=7.5),strip.text=element_text(size=8),axis.ticks=element_blank(),axis.line=element_blank(),panel.grid=element_blank(),panel.spacing=unit(2,'mm'),legend.position='bottom',legend.text=element_text(size=7.5),legend.title=element_text(size=8))+guides(fill=guide_colourbar(barwidth=unit(50,'mm'),barheight=unit(3,'mm')))
  save_panel(p,paste0('Scale_Matched_Ranks_',flow),185)
 }
}
writeLines('Reused original generate_main_spatial_validation_figure.R setup and base_theme. Arial throughout; vector PDF/SVG and 300-dpi PNG, width 183 mm. Original manuscript figures untouched.',file.path(out,'style_provenance.txt'))
