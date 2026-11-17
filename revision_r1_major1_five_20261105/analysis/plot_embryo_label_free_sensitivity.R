#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(dplyr);library(ggplot2);library(readr)})
source("external/atlas-feature-selection-benchmark/analysis/R/plotting.R")
args <- commandArgs(trailingOnly=TRUE)
src <- args[[1]]; out <- args[[2]]; dir.create(out,recursive=TRUE,showWarnings=FALSE)
settings <- read_tsv(file.path(src,"representative_settings.tsv"),show_col_types=FALSE,col_types=cols(n_features=col_character()))
metrics <- read_tsv(file.path(src,"scaled_setting_metrics.tsv"),show_col_types=FALSE,col_types=cols(n_features=col_character()))
method_labels <- c("scanpy_seurat_v3"="Scanpy Seurat v3","scanpy_seurat_v3_batch"="Scanpy Seurat v3 (batch)","scanpy_cell_ranger"="Scanpy Cell Ranger","scanpy_cell_ranger_batch"="Scanpy Cell Ranger (batch)","scanpy_seurat"="Scanpy Seurat","scanpy_seurat_batch"="Scanpy Seurat (batch)","scanpy_pearson"="Scanpy Pearson","scanpy_pearson_batch"="Scanpy Pearson (batch)","seurat_vst"="Seurat VST","seurat_mvp"="Seurat MVP","seurat_disp"="Seurat Dispersion","seurat_sct"="Seurat sctransform","singleCellHaystack"="singleCellHaystack","statistic_mean"="Mean expression","statistic_variance"="Expression variance","morans_i"="Moran's I","spatialde"="SpatialDE","sparkx"="SPARK-X","nnsvg"="nnSVG","somde"="SOMDE","Brennecke"="Brennecke","scPNMF"="scPNMF","nbumi"="NBumi","dubstepr"="DUBStepR","osca"="OSCA","scry"="scry","anticor"="Anticor","hotspot"="Hotspot","triku"="triku","all_features"="All features","random"="Random (N=500)")
settings <- settings |> filter(dataset %in% c("E8p5Embryo","E9p5Embryo")) |> filter(MethodGroup %in% c("expression_driven","spatially_informed") | fs_method %in% c("all_features","random")) |> select(dataset,fs_method,n_features,MethodGroup)
d <- metrics |> inner_join(settings,by=c("dataset","fs_method","n_features")) |> filter(metric_name %in% c("iLISI","CHAOS","PAS","silhouette"),integration_method %in% c("scvi","cellcharter")) |> mutate(
 Dataset=factor(dataset,levels=c("E8p5Embryo","E9p5Embryo"),labels=c("E8.5 embryo","E9.5 embryo")),
 Integrator=factor(integration_method,levels=c("scvi","cellcharter"),labels=c("scVI","CellCharter")),
 Metric=factor(metric_name,levels=c("iLISI","CHAOS","PAS","silhouette"),labels=c("iLISI","CHAOS","PAS","Silhouette")),
 Method=recode(fs_method,!!!method_labels), Type=if_else(fs_method %in% c("all_features","random"),"Reference","Competitive method"))
stopifnot(nrow(d)==496, n_distinct(d$fs_method)==31, !anyDuplicated(d[c("dataset","integration_method","fs_method","metric_name")]))
method_order <- d |> distinct(fs_method,Method,Type) |> arrange(Type,Method) |> pull(Method)
d <- d |> mutate(Method=factor(Method,levels=rev(method_order)))
p <- ggplot(d,aes(x=Metric,y=Method,fill=ScaledValue)) +
 geom_tile(colour="white",linewidth=0.15,width=.91,height=.91) +
 facet_grid(.~ Dataset + Integrator,scales="free_x",space="free_x") +
 scale_fill_viridis_c(option="magma",direction=1,limits=c(0,1),breaks=c(0,.5,1),na.value="grey85",name="Scaled value (oriented)",guide=guide_colourbar(barwidth=grid::unit(4,"cm"),title.position="top")) +
 labs(x=NULL,y=NULL) + theme_features_pub() +
 theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text=element_text(colour="black"),axis.text.x=element_text(angle=90,hjust=1,vjust=.5,size=6),axis.text.y=element_text(size=6),panel.grid=element_blank(),strip.background=element_rect(fill="black",colour="black",linewidth=.3),strip.text=element_text(colour="white",face="bold",size=7),panel.spacing.x=grid::unit(1.5,"mm"),legend.position="bottom",legend.title=element_text(size=7),legend.text=element_text(size=6),plot.margin=margin(4,4,4,4))
saveRDS(p,file.path(out,"embryo_label_free_sensitivity.rds"))
save_figure_files(p,file.path(out,"embryo_label_free_sensitivity"),width=8.2,height=8.0)
ggsave(file.path(out,"embryo_label_free_sensitivity.svg"),p,width=8.2,height=8.0,device=svglite::svglite,bg="white")
ggsave(file.path(out,"embryo_label_free_sensitivity.tiff"),p,width=8.2,height=8.0,dpi=600,device=ragg::agg_tiff,compression="lzw",bg="white")
write_tsv(d |> select(dataset,fs_method,n_features,MethodGroup,integration_method,metric_name,ValueMeanRaw,ValueSDRaw,OrientedMean,ScaledValue,NSeeds,EffectiveFeaturesMin,EffectiveFeaturesMax),file.path(out,"embryo_label_free_sensitivity_source.tsv"),na="NA")
cat("Wrote",nrow(d),"observations\n")
