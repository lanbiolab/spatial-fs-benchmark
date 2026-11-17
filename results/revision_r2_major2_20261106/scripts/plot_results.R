original <- readLines('scripts/generate_main_spatial_validation_figure.R')
setup_end <- grep('^semi <-', original)[1]-1L
eval(parse(text=paste(original[seq_len(setup_end)],collapse='\n')))
root <- 'results/revision_r2_major2_20261106'
output_dir <- file.path(root,'figures');dir.create(output_dir,showWarnings=FALSE)
d <- read_tsv(file.path(root,'analysis/individual_summary.tsv'),show_col_types=FALSE)
order <- c('seurat_vst','scanpy_seurat_v3_batch','triku','scanpy_cell_ranger_batch','singleCellHaystack','morans_i','sparkx','nnsvg','spatialde','somde','random','all_features','wilcoxon','scsegindex')
d <- d |> mutate(Method=factor(method_label(.data$method),levels=rev(method_label(order))))
draw <- function(x,xlabel,cols=2) {
 means <- x |> group_by(.data$Method,.data$Panel) |> summarise(Mean=mean(.data$Value),.groups='drop')
 bands <- x |> filter(.data$role=='reference') |> distinct(.data$Method,.data$Panel) |> mutate(y=as.numeric(.data$Method))
 ggplot(x,aes(x=.data$Value,y=.data$Method))+
 geom_rect(data=bands,aes(xmin=-Inf,xmax=Inf,ymin=.data$y-.47,ymax=.data$y+.47),inherit.aes=FALSE,fill='grey70',alpha=.28,colour=NA)+
 geom_point(aes(colour=.data$dataset),size=1.15,alpha=.65,position=position_jitter(height=.1,width=0,seed=17))+
 geom_point(data=means,aes(x=.data$Mean,y=.data$Method),inherit.aes=FALSE,shape=23,size=1.8,stroke=.3,fill='grey20',colour='white')+
 facet_wrap(~Panel,ncol=cols,scales='free_x')+scale_y_discrete(drop=FALSE)+
 scale_colour_manual(values=c(DLPFC='#377eb8',STOmics0212='#e41a1c',STOmics0218='#4daf4a',STOmics0224='#984ea3',E9p5Embryo='#ff7f00'),guide='none')+
 labs(x=xlabel,y=NULL)+base_theme+
 theme(axis.text.y=element_text(size=5.5),axis.text.x=element_text(size=5.5),axis.line.y=element_blank(),axis.ticks.y=element_blank(),strip.text=element_text(size=6,colour='white',face='bold'),panel.spacing=unit(.35,'cm'))
}
save <- function(p,name,h) {
 ggsave(file.path(output_dir,paste0(name,'.pdf')),p,width=183,height=h,units='mm',device=grDevices::cairo_pdf,family='Arial',bg='white')
 ggsave(file.path(output_dir,paste0(name,'.png')),p,width=183,height=h,units='mm',dpi=600,bg='white')
}
a <- d |> filter(.data$dataset %in% c('DLPFC','STOmics0212')) |> mutate(Panel=if_else(.data$dataset=='DLPFC','DLPFC (2 donors)','Skin / STOmics0212 (7 donors)'),Value=.data$heldout_ari)
save(draw(a,'Held-out individual ARI'), 'R2_M2_individual_ARI',106)
write_tsv(a,file.path(output_dir,'R2_M2_individual_ARI_source.tsv'))
b <- d |> filter(is.finite(.data$heldout_moran_percentile)) |> mutate(Panel=factor(.data$dataset,levels=c('DLPFC','STOmics0212','STOmics0218','STOmics0224','E9p5Embryo')),Value=.data$heldout_moran_percentile)
save(draw(b,'Mean held-out Moran percentile',3)+geom_vline(xintercept=.5,linetype='dashed',colour='grey60',linewidth=.28),'R2_M2_individual_spatial',166)
write_tsv(b,file.path(output_dir,'R2_M2_individual_spatial_source.tsv'))
