suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(ggforce);library(colorspace)})
out='results/revision_four_item_20261109'
save=function(p,n,w,h){set.seed(20261109);ggsave(file.path(out,paste0(n,'.pdf')),p,width=w/72,height=h/72,units='in',device=grDevices::cairo_pdf,family='Arial',bg='white')}
check=function(a,b){stopifnot(identical(a$data,b$data),identical(a$mapping,b$mapping),length(a$layers)==length(b$layers));for(j in seq_along(a$layers))stopifnot(identical(a$layers[[j]]$layer_data(a$data),b$layers[[j]]$layer_data(b$data)),identical(a$layers[[j]]$mapping,b$layers[[j]]$mapping))}
p=readRDS('results/revision_all_figure_readability_20261109/Figure2.rds'); old=p
p[[1]][[2]]=p[[1]][[2]]+theme(strip.text.x=element_text(size=6.8,angle=90),strip.text.y=element_text(size=6.8))
p[[2]][[2]]=p[[2]][[2]]+theme(axis.text.x=element_text(size=6.2))
for(j in seq_along(p[[2]][[2]]$layers))if(inherits(p[[2]][[2]]$layers[[j]]$geom,'GeomText'))p[[2]][[2]]$layers[[j]]$aes_params$size=6.2/ggplot2::.pt
for(i in 1:2)for(j in seq_len(length(p[[i]])))check(old[[i]][[j]],p[[i]][[j]])
save(p,'Figure2',540,565)
for(n in c('b','c','d')){
 p=readRDS(paste0('results/revision_all_figure_readability_20261109/Figure4',n,'.rds'));old=p
 cat(n,'coord',class(p$coordinates),'aspect',p$theme$aspect.ratio,'\n')
 for(k in c('axis.text.x','axis.text.y','strip.text','legend.text','legend.title')){
  e=ggplot2::calc_element(k,p$theme)
  if(!inherits(e,'element_blank')){e$size=if(k=='legend.title')6.5 else 6.2;p=p+do.call(theme,setNames(list(e),k))}
 }
 p=p+theme(plot.tag.position=c(.015,.98),plot.tag=element_text(size=9,face='bold',vjust=1),plot.margin=margin(4,3,3,3))
 check(old,p);save(p,paste0('Figure4',n),c(b=252,c=126,d=162)[n],230)
}
writeLines('Original data, mappings, layer data and ordering unchanged for all edited plot objects.',file.path(out,'plot_verification.txt'))
