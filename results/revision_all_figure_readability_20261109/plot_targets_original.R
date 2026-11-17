original <- readLines('scripts/generate_main_spatial_validation_figure.R')
setup_end <- grep('^semi <-',original)[1]-1L
eval(parse(text=paste(original[seq_len(setup_end)],collapse='\n')))
root <- 'results/revision_r2_major4_20261107'
output_dir <- file.path(root,'figures');dir.create(output_dir,showWarnings=FALSE)
d <- read_tsv(file.path(root,'donor_summary.tsv'),show_col_types=FALSE) |>
 mutate(strategy=factor(strategy,levels=c('HVG','SVG','Hybrid')),condition=factor(condition,levels=c('Unconstrained','Retain targets')))
draw <- function(metric,label,limits=NULL) {
 p <- ggplot(d,aes(x=condition,y=.data[[metric]],group=heldout_donor,colour=heldout_donor))+
  geom_line(linewidth=.38,alpha=.75)+geom_point(aes(shape=heldout_donor),size=1.8)+
  facet_wrap(~strategy,nrow=1)+scale_colour_manual(values=c(Br5292='#377eb8',Br8100='#e41a1c'),name='Held-out donor')+
  scale_shape_manual(values=c(Br5292=16,Br8100=17),name='Held-out donor')+
  labs(x=NULL,y=label)+base_theme+
  theme(axis.text.x=element_text(size=6.5),axis.text.y=element_text(size=6.5),axis.title.y=element_text(size=7),strip.text=element_text(size=7),panel.spacing=unit(.28,'cm'),legend.position='bottom',legend.text=element_text(size=7),legend.title=element_text(size=7))
 if(!is.null(limits))p <- p+coord_cartesian(ylim=limits)
 p
}
p <- (draw('target_retention','Target retention',c(0,1.04))/draw('heldout_ari','Held-out ARI')/draw('heldout_nmi','Held-out NMI')/draw('heldout_moran_percentile','Held-out Moran percentile'))+
 plot_layout(guides='collect')+plot_annotation(tag_levels='a')&theme(legend.position='bottom')
ggsave(file.path(output_dir,'R2_M4_target_retention.pdf'),p,width=183,height=166,units='mm',device=grDevices::cairo_pdf,family='Arial',bg='white')
ggsave(file.path(output_dir,'R2_M4_target_retention.png'),p,width=183,height=166,units='mm',dpi=600,bg='white')
write_tsv(d,file.path(output_dir,'R2_M4_target_retention_source.tsv'))
