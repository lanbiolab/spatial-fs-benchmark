original <- readLines('scripts/generate_main_spatial_validation_figure.R')
eval(parse(text=paste(original[seq_len(grep('^semi <-',original)[1]-1L)],collapse='\n')))
root <- 'results/revision_r3_comment5_20261107'
output_dir <- file.path(root,'figures');dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
a <- file.path(root,'analysis')
d <- read_tsv(file.path(a,'03_competitive_endpoint_ranks.tsv'),show_col_types=FALSE)
order <- sort(unique(d$fs_method))
labs <- method_label(order)
labs[order %in% c('dubstepr','nnsvg','seurat_mvp')] <- paste0(labs[order %in% c('dubstepr','nnsvg','seurat_mvp')],' *')
platforms <- c('visium','st','stereo_seq','openst');pn <- c('Visium','ST','Stereo-seq','Open-ST')
d <- d |> mutate(Method=factor(fs_method,levels=rev(order),labels=rev(labs)),Platform=factor(platform_id,levels=platforms,labels=pn),Metric=factor(metric_name,levels=c('silhouette','CHAOS','PAS','ari','nmi'),labels=c('Silhouette (higher)','CHAOS (lower)','PAS (lower)','ARI (higher)','NMI (higher)')))
save <- function(p,name,h) {
 ggsave(file.path(output_dir,paste0(name,'.pdf')),p,width=183,height=h,units='mm',device=grDevices::cairo_pdf,family='Arial',bg='white')
 ggsave(file.path(output_dir,paste0(name,'.png')),p,width=183,height=h,units='mm',dpi=300,device=ragg::agg_png,bg='white')
 ggsave(file.path(output_dir,paste0(name,'.tiff')),p,width=183,height=h,units='mm',dpi=600,device=ragg::agg_tiff,compression='lzw',bg='white')
 ggsave(file.path(output_dir,paste0(name,'.svg')),p,width=183,height=h,units='mm',device=svglite::svglite,bg='white')
}
for (workflow in c('scvi','cellcharter','graphst')) {
 x <- d |> filter(integration_method==workflow,metric_name %in% c('silhouette','CHAOS','PAS'))
 p <- ggplot(x,aes(Platform,Method,fill=rank)) + geom_tile(colour='white',linewidth=.15) +
  geom_point(data=x |> filter(rank<=3),shape=21,fill='white',colour='black',stroke=.25,size=.9)+
  facet_grid(.~Metric)+scale_fill_gradientn(colours=c('#2166ac','#67a9cf','#d1e5f0','#f7f7f7'),limits=c(1,29),breaks=c(1,10,20,29),name='Rank (1 = best)')+
  labs(x=NULL,y=NULL)+base_theme+theme(axis.text.y=element_text(size=6),axis.text.x=element_text(size=6,angle=40,hjust=1),axis.ticks=element_blank(),axis.line=element_blank(),panel.grid=element_blank(),strip.text=element_text(size=7),legend.position='bottom',legend.text=element_text(size=6),legend.title=element_text(size=7),legend.key.height=unit(2,'mm'),legend.key.width=unit(4,'mm'),panel.spacing=unit(.15,'cm'))
 save(p,paste0('R3C5_common_ranks_',workflow),156)
}
s <- read_tsv(file.path(a,'02_metric_means_seed_sd.tsv'),show_col_types=FALSE) |>
 filter(platform_id %in% c('st','stereo_seq'),metric_name=='ari')
refs <- c('TFs','random','all_features','scsegindex','wilcoxon');all_order <- c(order,refs)
s <- s |> mutate(Method=factor(fs_method,levels=rev(all_order),labels=rev(method_label(all_order))),Platform=factor(platform_id,levels=c('st','stereo_seq'),labels=c('ST: atlas-registered reference','Stereo-seq: source annotation')),Consumer=factor(integration_method,levels=c('scvi','cellcharter','graphst'),labels=c('scVI','CellCharter','GraphST')),Low=pmax(-1,mean-sd),High=pmin(1,mean+sd))
bands <- s |> filter(role=='reference') |> distinct(Method,Platform,Consumer) |> mutate(y=as.numeric(Method))
p <- ggplot(s,aes(mean,Method))+
 geom_rect(data=bands,aes(xmin=-Inf,xmax=Inf,ymin=y-.48,ymax=y+.48),inherit.aes=FALSE,fill='grey70',alpha=.28)+
 geom_linerange(aes(xmin=Low,xmax=High),linewidth=.35,colour='grey35')+
 geom_point(aes(colour=role),size=1.25)+facet_grid(Platform~Consumer)+
 scale_y_discrete(drop=FALSE)+scale_colour_manual(values=c(competitive='#377eb8',reference='#7A5AA6'),guide='none')+
 scale_x_continuous(limits=c(-.1,.85),breaks=c(0,.25,.5,.75))+
 labs(x='ARI against supplied reference (mean +/- seed SD)',y=NULL)+base_theme+
 theme(axis.text.y=element_text(size=5.1),axis.text.x=element_text(size=6),strip.text=element_text(size=6.3),panel.spacing=unit(.15,'cm'),axis.ticks.y=element_blank(),axis.line.y=element_blank())
save(p,'R3C5_annotation_ARI',215)
r <- read_tsv(file.path(a,'06_cross_platform_rank_correlations.tsv'),show_col_types=FALSE)
pair_order <- c('Visium / ST','Visium / Stereo-seq','Visium / Open-ST','ST / Stereo-seq','ST / Open-ST','Stereo-seq / Open-ST')
r <- r |> mutate(Pair=factor(paste(pn[match(platform1,platforms)],pn[match(platform2,platforms)],sep=' / '),levels=pair_order),Condition=factor(paste(workflow,metric,sep=': '),levels=rev(as.vector(t(outer(c('scvi','cellcharter','graphst'),c('silhouette','CHAOS','PAS'),paste,sep=': '))))),Subset=factor(subset,levels=c('all_29','exact_2000_all_four'),labels=c('All 29 competitive selectors','26 selectors with exactly 2,000 genes')))
r$Condition <- factor(gsub('graphst:', 'GraphST:', gsub('cellcharter:', 'CellCharter:', gsub('scvi:', 'scVI:', as.character(r$Condition)))),levels=gsub('graphst:', 'GraphST:', gsub('cellcharter:', 'CellCharter:', gsub('scvi:', 'scVI:', levels(r$Condition)))))
p <- ggplot(r,aes(Pair,Condition,fill=rho))+geom_tile(colour='white',linewidth=.3)+geom_text(aes(label=sprintf('%.2f',rho)),size=2,colour='black')+
 facet_grid(.~Subset)+scale_fill_gradient2(low='#b2182b',mid='white',high='#2166ac',limits=c(-1,1),midpoint=0,name='Spearman rho')+
 labs(x=NULL,y=NULL)+guides(fill=guide_colourbar(barwidth=unit(42,'mm'),barheight=unit(2.5,'mm')))+base_theme+theme(axis.text.y=element_text(size=6),axis.text.x=element_text(size=5.8,angle=55,hjust=1),strip.text=element_text(size=6.7),axis.ticks=element_blank(),axis.line=element_blank(),legend.position='bottom',legend.key.height=unit(2,'mm'),legend.text=element_text(size=6),legend.title=element_text(size=7),panel.spacing=unit(.25,'cm'))
save(p,'R3C5_platform_correlations',104)
writeLines('Original generate_main_spatial_validation_figure.R setup and theme_features_pub() reused without edits. New supplemental panels only; 183 mm width, Arial family, PDF/SVG/600-dpi TIFF/300-dpi PNG.',file.path(output_dir,'style_provenance.txt'))
