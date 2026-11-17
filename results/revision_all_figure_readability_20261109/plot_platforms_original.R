original <- readLines("scripts/generate_main_spatial_validation_figure.R")
eval(parse(text=paste(original[seq_len(grep("^semi <-",original)[1]-1L)],collapse="\n")))
root <- "results/revision_r2_major1_20261106"
analysis <- file.path(root,"analysis")
dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
tab <- read_tsv(file.path(analysis,"platform_common_component_method_summary.tsv"),show_col_types=FALSE)
wide <- tab |>
  mutate(Platform=if_else(grepl("Visium",platform),"Visium","SlideSeq")) |>
  select(integration_method,Component,fs_method,Platform,PlatformRank) |>
  pivot_wider(names_from=Platform,values_from=PlatformRank) |>
  mutate(Integrator=factor(integration_method,levels=c("scvi","cellcharter"),labels=c("scVI","CellCharter")),
         Task=factor(Component,levels=c("BatchMixing","LabelIndependentClustering"),labels=c("Batch mixing","Label-independent clustering")),
         Class=if_else(fs_method %in% c("morans_i","sparkx","nnsvg","spatialde","somde"),"Spatially informed","Expression driven"),
         Label=if_else(fs_method %in% c("morans_i","nnsvg","somde","triku","seurat_vst"),method_label(fs_method),""))
p <- ggplot(wide,aes(Visium,SlideSeq)) +
  geom_abline(slope=1,intercept=0,colour="grey75",linetype="dashed",linewidth=.35) +
  geom_point(aes(colour=Class),size=1.7,alpha=.8) +
  facet_grid(Task~Integrator) +
  scale_colour_manual(values=c("Expression driven"="#377eb8","Spatially informed"="#4daf4a"),name=NULL) +
  scale_x_continuous(breaks=c(1,5,10,15,20,25,29),limits=c(0,30)) +
  scale_y_continuous(breaks=c(1,5,10,15,20,25,29),limits=c(0,30)) +
  labs(x="Visium-only within-platform rank (five datasets)",y="Slide-seqV2 within-platform rank (two stages, one study)") +
  base_theme + theme(axis.text=element_text(size=7),axis.title=element_text(size=8),
                     strip.text=element_text(size=8),legend.text=element_text(size=8),legend.position="bottom")
if(requireNamespace("ggrepel",quietly=TRUE)){
  p <- p + ggrepel::geom_text_repel(aes(label=Label),size=2.4,min.segment.length=0,seed=20261106,max.overlaps=Inf,box.padding=.2)
} else p <- p+geom_text(aes(label=Label),size=2.2,check_overlap=TRUE,vjust=-.5)

r <- read_tsv(file.path(analysis,"conditional_physical_radius_neighborhood_audit.tsv"),show_col_types=FALSE) |>
  filter(radius_um_assumed==150) |>
  mutate(Dataset=recode(dataset,"MouseBrainSerialSections"="Mouse Brain","E9p5Embryo"="E9.5 (conditional)",
                        "STOmics0212"="STOmics-0212","STOmics0218"="STOmics-0218","STOmics0224"="STOmics-0224"),
         Dataset=factor(Dataset,levels=c("DLPFC","Mouse Brain","STOmics-0212","STOmics-0218","STOmics-0224","E9.5 (conditional)"))) |>
  pivot_longer(c(knn8_support_median_um_assumed,degree_median),names_to="Endpoint",values_to="Value") |>
  mutate(Endpoint=factor(Endpoint,levels=c("knn8_support_median_um_assumed","degree_median"),
                         labels=c("Median distance to eighth neighbour (conditional um)","Median neighbour count within 150 um (conditional)")))
q <- ggplot(r,aes(Dataset,Value,colour=Dataset)) +
  geom_point(position=position_jitter(width=.10,height=0,seed=20261106),size=1.6,alpha=.75) +
  facet_wrap(~Endpoint,scales="free_y",nrow=1) +
  scale_colour_manual(values=c("#377eb8","#4daf4a","#e41a1c","#984ea3","#ff7f00","#a65628"),guide="none") +
  labs(x=NULL,y=NULL) + base_theme +
  theme(axis.text.x=element_text(size=7,angle=40,hjust=1),axis.text.y=element_text(size=7),
        strip.text=element_text(size=7),plot.margin=margin(3,3,3,3,"mm"))
save_candidates <- function(plot,name,w,h){
  ggsave(file.path(output_dir,paste0(name,".pdf")),plot,width=w,height=h,device=grDevices::cairo_pdf,family="Arial",bg="white")
  ggsave(file.path(output_dir,paste0(name,".png")),plot,width=w,height=h,dpi=600,device=ragg::agg_png,bg="white")
  ggsave(file.path(output_dir,paste0(name,".svg")),plot,width=w,height=h,device=svglite::svglite,bg="white")
}
save_candidates(p,"platform_common_component_ranks",7.8,5.8)
save_candidates(q,"neighborhood_support_diagnostic",7.8,3.8)
write_tsv(wide,file.path(output_dir,"platform_rank_plot_source.tsv"))
write_tsv(r,file.path(output_dir,"neighborhood_plot_source.tsv"))
message("New R2-M1 diagnostic candidates exported; original figures were not modified.")
