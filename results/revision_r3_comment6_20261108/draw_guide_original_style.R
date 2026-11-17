suppressPackageStartupMessages({
 library(ggplot2); library(dplyr); library(tidyr); library(patchwork); library(cowplot)
})
source("external/atlas-feature-selection-benchmark/analysis/R/plotting.R")
register_arial_pdf_font()
out <- "results/revision_r3_comment6_20261108/figures"
dir.create(out,recursive=TRUE,showWarnings=FALSE)
tasks <- c("Spatial-pattern discovery","Expression-only clustering","Spatially informed clustering","Multi-slice integration","Spatial alignment")
methods <- c("Moran's I","nnSVG","Seurat-VST","scanpy-SeuratV3 (batch)","triku","HVG-SVG hybrid")
cols <- setNames(c("#8073ac","#377eb8","#377eb8","#e41a1c","#4daf4a"),tasks)
chosen <- list(c(1,2),c(3,4,5),c(1,3,4,6),c(1,3,4),c(1,3))
d <- bind_rows(lapply(seq_along(tasks),function(i) data.frame(Task=tasks[i],Method=methods[chosen[[i]]],Role=if(i<=2)"Starting shortlist" else "Compare in workflow")))
d$Task <- factor(d$Task,levels=rev(tasks));d$Method <- factor(d$Method,levels=methods)
d$Panel <- "Candidate selectors"
theme_integration <- theme_features_pub() +
 theme(axis.title=element_blank(),
       axis.text.x=element_text(size=5.0,angle=90,hjust=1,vjust=.5,face="plain",colour="black"),
       axis.text.y=element_text(size=5.1,face="plain",colour="black"),
       panel.grid=element_blank(),
       strip.background=element_rect(fill="black",colour="black"),
       strip.text=element_text(colour="white",face="bold"),
       legend.position="bottom",
       legend.title=element_text(size=6,face="bold"),
       legend.text=element_text(size=5),
       plot.margin=margin(2,2,2,2))
p <- ggplot(d,aes(Method,Task,colour=Task,shape=Role))+
 geom_point(size=2.05,stroke=.55)+
 facet_grid(.~Panel)+
 scale_colour_manual(values=cols,guide="none")+
 scale_shape_manual(values=c("Starting shortlist"=15,"Compare in workflow"=0),
                    breaks=c("Starting shortlist","Compare in workflow"))+
 scale_x_discrete(drop=FALSE,expand=expansion(add=.5))+
 scale_y_discrete(drop=FALSE,expand=expansion(add=.5))+
 guides(shape=guide_legend(title="Candidate role",nrow=1,override.aes=list(colour="black",size=2.05)))+
 theme_integration
crit <- data.frame(Task=factor(tasks,levels=rev(tasks)),
 Text=c("External or held-out references\nSpatial signal and biological relevance",
 "Annotation agreement and stability\nWithin the same embedding workflow",
 "Annotation agreement and spatial coherence\nEvaluate the two endpoints separately",
 "Batch mixing and biological conservation\nUse biologically meaningful batches",
 "Spatial correspondence and anatomical agreement\nValidate on relevant tissue pairs"),
 Panel="Task-specific decision criteria")
q <- ggplot(crit,aes(0,Task,label=Text))+
 geom_text(hjust=0,size=5.3/ggplot2::.pt,family="Arial",lineheight=1.15)+
 facet_grid(.~Panel)+
 scale_x_continuous(limits=c(0,1),expand=expansion(mult=c(.035,.015)))+
 scale_y_discrete(drop=FALSE,expand=expansion(add=.5))+
 theme_integration+
 theme(axis.text.x=element_blank(),axis.text.y=element_blank(),
       axis.ticks.x=element_blank(),axis.ticks.y=element_blank(),legend.position="none")
main <- (p+q+plot_layout(widths=c(1.15,1)))+
 plot_annotation(tag_levels=list(c("e","")))&
 theme(plot.tag=element_text(size=9,face="bold",family="Arial"),plot.tag.position=c(0,1))
foot <- ggdraw()+
 draw_label("Required genes: reserve within budget if they must shape the representation; retain the full matrix for biological testing.",
            x=.155,y=.81,hjust=0,fontfamily="Arial",size=5.3)+
 draw_label("Match actual gene counts and workflow; compare stability and compute cost. Blank cells do not imply an unsuitable method.",
            x=.155,y=.35,hjust=0,fontfamily="Arial",size=5.3)
fig <- plot_grid(main,foot,ncol=1,rel_heights=c(1,.105))
save_figure_files(fig,file.path(out,"Figure6e_task_guide"),width=8.2,height=3.1)
svglite::svglite(file.path(out,"Figure6e_task_guide.svg"),width=8.2,height=3.1,bg="white");print(fig);dev.off()
ragg::agg_tiff(file.path(out,"Figure6e_task_guide.tiff"),width=8.2,height=3.1,units="in",res=600,compression="lzw",background="white");print(fig);dev.off()
write.table(d,file.path(out,"guidance_matrix.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
saveRDS(fig,file.path(out,"Figure6e_task_guide.rds"))
writeLines(capture.output(warnings()),file.path(out,"plot_warnings.txt"))

