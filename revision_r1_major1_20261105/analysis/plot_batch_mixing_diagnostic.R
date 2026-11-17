#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(dplyr);library(ggplot2);library(readr);library(tidyr)})
source("external/atlas-feature-selection-benchmark/analysis/R/plotting.R")
args <- commandArgs(trailingOnly = TRUE)
d <- read_tsv(args[[1]], show_col_types=FALSE)
out <- args[[2]]; dir.create(out,recursive=TRUE,showWarnings=FALSE)
d <- d |> filter(MethodGroup %in% c("expression_driven","spatially_informed") | fs_method %in% c("random","all_features")) |>
 mutate(Dataset = factor(dataset, levels=c("DLPFC","STOmics0212","STOmics0218"), labels=c("DLPFC","STOmics-0212","STOmics-0218")),
 Integrator=factor(integration_method,levels=c("scvi","cellcharter"),labels=c("scVI","CellCharter")),
 Group=case_when(fs_method=="random" ~ "Random (500 genes)",fs_method=="all_features" ~ "All features",TRUE ~ "Competitive methods (2,000 requested)"))
comp <- d |> filter(Group=="Competitive methods (2,000 requested)")
ctrl <- d |> filter(Group!="Competitive methods (2,000 requested)")
med <- comp |> group_by(Dataset,Integrator) |> summarise(BatchMixingMean=median(BatchMixingMean,na.rm=TRUE),.groups="drop")
p <- ggplot(comp,aes(x=Dataset,y=BatchMixingMean)) +
 geom_point(position=position_jitter(width=.14,height=0,seed=4),shape=16,size=1.1,alpha=.45,colour="grey48") +
 geom_point(data=med,shape=95,size=6,colour="black") +
 geom_point(data=ctrl,aes(colour=Group,shape=Group),position=position_dodge(width=.28),size=2.7,stroke=.6) +
 facet_wrap(~Integrator,nrow=1) +
 scale_colour_manual(values=c("Random (500 genes)"="#f781bf","All features"="#377eb8"),name=NULL) +
 scale_shape_manual(values=c("Random (500 genes)"=18,"All features"=17),name=NULL) +
 labs(x=NULL,y="Batch-mixing component (scaled)",caption="Gray points: competitive methods; black bar: competitive median") +
 theme_features_pub() + theme(panel.grid.major.x=element_blank(),panel.grid.minor=element_blank(),panel.border=element_rect(fill=NA,colour="black",linewidth=.3),strip.background=element_rect(fill="black",colour="black"),strip.text=element_text(colour="white",face="bold"),axis.text.x=element_text(angle=30,hjust=1),legend.position="bottom",plot.caption=element_text(size=7,hjust=0))
saveRDS(p,file.path(out,"batch_mixing_primary.rds"))
save_figure_files(p,file.path(out,"batch_mixing_primary"),width=7.2,height=3.8)
write_tsv(d,file.path(out,"batch_mixing_primary_source.tsv"),na="NA")

ggsave(file.path(out,"batch_mixing_primary.svg"),p,width=7.2,height=3.8,device=svglite::svglite,bg="white")
ggsave(file.path(out,"batch_mixing_primary.tiff"),p,width=7.2,height=3.8,dpi=600,device=ragg::agg_tiff,compression="lzw",bg="white")
