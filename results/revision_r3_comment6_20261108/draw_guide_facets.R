suppressPackageStartupMessages({
 library(ggplot2); library(grid)
})
source("external/atlas-feature-selection-benchmark/analysis/R/plotting.R")
register_arial_pdf_font()
out <- "results/revision_r3_comment6_20261108/figures"
dir.create(out,recursive=TRUE,showWarnings=FALSE)
sty <- theme_features_pub()
W<-590;H<-154
g<-list()
add<-function(x){g[[length(g)+1]]<<-x}
txt<-function(x,y,s,size=5.1,bold=FALSE,col="black",just="left",lineheight=1.1){
 add(textGrob(s,x=unit(x,"bigpts"),y=unit(H-y,"bigpts"),just=c(just,"centre"),
 gp=gpar(fontfamily=sty$text$family,fontsize=size,fontface=if(bold)"bold" else "plain",col=col,lineheight=lineheight)))
}
box<-function(x,y,w,h,fill=NA,col="black",lwd=.65){
 add(rectGrob(x=unit(x,"bigpts"),y=unit(H-y,"bigpts"),width=unit(w,"bigpts"),height=unit(h,"bigpts"),just=c("left","top"),gp=gpar(fill=fill,col=col,lwd=lwd)))
}
ln<-function(x,y,x2,y2,col="black",lwd=.65){
 add(segmentsGrob(x0=unit(x,"bigpts"),y0=unit(H-y,"bigpts"),x1=unit(x2,"bigpts"),y1=unit(H-y2,"bigpts"),gp=gpar(col=col,lwd=lwd)))
}
sq<-function(x,y,col,filled=TRUE,size=4.6){
 box(x-size/2,y-size/2,size,size,fill=if(filled)col else "white",col=col,lwd=.8)
}
tasks<-c("Spatial-pattern\ndiscovery","Expression-only\nclustering","Spatially informed\nclustering","Multi-slice\nintegration","Spatial\nalignment")
methods<-c("triku","Seurat-VST","scanpy-SeuratV3 (batch)","nnSVG","Moran's I","HVG-SVG hybrid")
cols<-c("#8073ac","#377eb8","#377eb8","#e41a1c","#4daf4a")
chosen<-list(c(4,5),c(1,2,3),c(2,3,5,6),c(2,3,5),c(2,5))
checks<-c("External / held-out references\nSpatial and biological relevance",
"Annotation agreement\nRepresentation stability",
"Annotation agreement\nSpatial coherence (separately)",
"Batch mixing\nBiological conservation",
"Spatial correspondence\nAnatomical agreement")
txt(3,7,"e",9,TRUE)
xleft<-98; gap<-2; fw<-(590-xleft-2-4*gap)/5
header_y<-3;header_h<-18
row_y<-seq(31,86,length.out=6)
for(i in seq_along(methods))txt(xleft-5,row_y[i],methods[i],5.1,just="right")
txt(xleft-5,110,"Select using",5.1,TRUE,just="right")
for(j in 1:5){
 x<-xleft+(j-1)*(fw+gap); cx<-x+fw/2
 box(x,header_y,fw,119,"white","black",.65)
 box(x,header_y,fw,header_h,"black","black",.65)
 txt(cx,12,tasks[j],sty$strip.text$size,TRUE,"white","centre",1.04)
 for(i in chosen[[j]])sq(cx,row_y[i],cols[j],j<=2)
 ln(x,96,x+fw,96,"black",.65)
 txt(cx,109.5,checks[j],5.1,FALSE,"black","centre",1.22)
}
txt(98,132.5,"Task",6,TRUE)
sq(158,132.5,cols[1],TRUE,4.4);txt(164,132.5,"Spatial-pattern discovery",5.1)
sq(248,132.5,cols[2],TRUE,4.4);txt(254,132.5,"Clustering",5.1)
sq(336,132.5,cols[4],TRUE,4.4);txt(342,132.5,"Integration",5.1)
sq(412,132.5,cols[5],TRUE,4.4);txt(418,132.5,"Alignment",5.1)
txt(98,144.5,"Candidate role",6,TRUE)
sq(158,144.5,"black",TRUE,4.4);txt(164,144.5,"Starting shortlist",5.1)
sq(248,144.5,"black",FALSE,4.4);txt(254,144.5,"Compare in workflow",5.1)
txt(336,144.5,"Blank: not shortlisted, not excluded",5.1)
fig<-do.call(grobTree,g)
save_figure_files(fig,file.path(out,"Figure6e_task_guide"),width=W/72,height=H/72)
svglite::svglite(file.path(out,"Figure6e_task_guide.svg"),width=W/72,height=H/72,bg="white");grid.draw(fig);dev.off()
ragg::agg_tiff(file.path(out,"Figure6e_task_guide.tiff"),width=W/72,height=H/72,units="in",res=600,compression="lzw",background="white");grid.draw(fig);dev.off()
saveRDS(fig,file.path(out,"Figure6e_task_guide.rds"))
writeLines(capture.output(warnings()),file.path(out,"plot_warnings.txt"))

