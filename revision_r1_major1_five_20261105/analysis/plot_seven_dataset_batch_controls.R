#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(dplyr);library(ggplot2);library(readr);library(tidyr)})
source("external/atlas-feature-selection-benchmark/analysis/R/plotting.R")
args <- commandArgs(trailingOnly=TRUE); src <- args[[1]]; out <- args[[2]]
dir.create(out,recursive=TRUE,showWarnings=FALSE)
settings <- read_tsv(file.path(src,"representative_settings.tsv"),show_col_types=FALSE,col_types=cols(n_features=col_character())) |>
 filter(MethodGroup %in% c("expression_driven","spatially_informed") | fs_method %in% c("random","all_features")) |>
 select(dataset,fs_method,n_features,MethodGroup)
tasks <- read_tsv(file.path(src,"representative_task_scores.tsv"),show_col_types=FALSE,col_types=cols(n_features=col_character())) |>
 semi_join(settings,by=c("dataset","fs_method","n_features")) |>
 transmute(dataset,fs_method,n_features,integration_method,MethodGroup,measure="Batch-mixing component",value=BatchMixingMean,n_seeds=BatchMixingNSeeds)
ilisi <- read_tsv(file.path(src,"scaled_setting_metrics.tsv"),show_col_types=FALSE,col_types=cols(n_features=col_character())) |>
 inner_join(settings,by=c("dataset","fs_method","n_features")) |>
 filter(metric_name=="iLISI") |>
 transmute(dataset,fs_method,n_features,integration_method,MethodGroup,measure="iLISI only",value=ScaledValue,n_seeds=NSeeds)
d <- bind_rows(tasks,ilisi) |> filter(integration_method %in% c("scvi","cellcharter")) |>
 mutate(Dataset=factor(dataset,levels=c("DLPFC","MouseBrainSerialSections","STOmics0212","STOmics0218","STOmics0224","E8p5Embryo","E9p5Embryo"),labels=c("DLPFC","Mouse Brain","STOmics-0212","STOmics-0218","STOmics-0224","E8.5*","E9.5*")),
 Integrator=factor(integration_method,levels=c("scvi","cellcharter"),labels=c("scVI","CellCharter")),
 Measure=factor(measure,levels=c("Batch-mixing component","iLISI only")),
 Group=case_when(fs_method=="random" ~ "Random (500 genes)",fs_method=="all_features" ~ "All features",TRUE ~ "Competitive methods"))
stopifnot(nrow(d)==868,all(d$n_seeds==3),!anyNA(d$value),!anyDuplicated(d[c("dataset","integration_method","measure","fs_method")]))
comp <- d |> filter(Group=="Competitive methods");ctrl <- d |> filter(Group!="Competitive methods")
stopifnot(all((comp |> count(dataset,integration_method,measure))$n==29),all((ctrl |> count(dataset,integration_method,measure))$n==2))
med <- comp |> group_by(Dataset,Integrator,Measure) |> summarise(value=median(value),.groups="drop")
p <- ggplot(comp,aes(x=Dataset,y=value)) +
 geom_point(position=position_jitter(width=.14,height=0,seed=4),shape=16,size=1.0,alpha=.42,colour="grey48") +
 geom_point(data=med,shape=95,size=5,colour="black") +
 geom_point(data=ctrl,aes(colour=Group,shape=Group),position=position_dodge(width=.32),size=2.35,stroke=.45) +
 facet_grid(Measure~Integrator) +
 scale_colour_manual(values=c("Random (500 genes)"="#f781bf","All features"="#377eb8"),name=NULL) +
 scale_shape_manual(values=c("Random (500 genes)"=18,"All features"=17),name=NULL) +
 scale_y_continuous(limits=c(-.02,1.02),breaks=c(0,.25,.5,.75,1)) +
 labs(x=NULL,y="Frozen-scaled score",caption="Gray points: 29 competitive methods; black bars: competitive medians") +
 theme_features_pub() + theme(panel.grid.major.x=element_blank(),panel.grid.minor=element_blank(),panel.border=element_rect(fill=NA,colour="black",linewidth=.3),strip.background=element_rect(fill="black",colour="black"),strip.text=element_text(colour="white",face="bold",size=7),axis.text.x=element_text(angle=40,hjust=1,size=6.5),axis.text.y=element_text(size=7),axis.title.y=element_text(size=8),legend.position="bottom",legend.text=element_text(size=7),plot.caption=element_text(size=6.5,hjust=0),panel.spacing=grid::unit(2,"mm"))
saveRDS(p,file.path(out,"batch_mixing_controls_seven.rds"))
save_figure_files(p,file.path(out,"batch_mixing_controls_seven"),width=8.2,height=6.0)
ggsave(file.path(out,"batch_mixing_controls_seven.svg"),p,width=8.2,height=6.0,device=svglite::svglite,bg="white")
ggsave(file.path(out,"batch_mixing_controls_seven.tiff"),p,width=8.2,height=6.0,dpi=600,device=ragg::agg_tiff,compression="lzw",bg="white")
write_tsv(d |> select(dataset,integration_method,measure,fs_method,n_features,MethodGroup,value,n_seeds),file.path(out,"batch_mixing_controls_seven_source.tsv"),na="NA")
summary <- ctrl |> rowwise() |> mutate(competitive_median=median(comp$value[comp$dataset==dataset & comp$integration_method==integration_method & comp$measure==measure]),rank_among_29_plus_control=rank(-c(value,comp$value[comp$dataset==dataset & comp$integration_method==integration_method & comp$measure==measure]),ties.method="average")[1]) |> ungroup() |>
 select(dataset,integration_method,measure,fs_method,n_features,value,competitive_median,rank_among_29_plus_control,n_seeds)
write_tsv(summary,file.path(out,"batch_mixing_controls_seven_summary.tsv"),na="NA")
cat("Completed",nrow(d),"metric observations;",nrow(summary),"reference comparisons\n")
