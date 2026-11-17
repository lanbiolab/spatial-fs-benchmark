suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(ggforce);library(colorspace)})
out='results/revision_layout_recheck_20261109'
save=function(p,n,w,h){set.seed(20261109);ggsave(file.path(out,paste0(n,'.pdf')),p,width=w/72,height=h/72,units='in',device=grDevices::cairo_pdf,family='Arial',bg='white')}
check=function(a,b){stopifnot(identical(a$data,b$data),identical(a$mapping,b$mapping),length(a$layers)==length(b$layers));for(j in seq_along(a$layers))stopifnot(identical(a$layers[[j]]$layer_data(a$data),b$layers[[j]]$layer_data(b$data)),identical(a$layers[[j]]$mapping,b$layers[[j]]$mapping))}
p=readRDS('results/revision_all_figure_readability_20261109/Figure2.rds');old=p
p[[1]][[2]]=p[[1]][[2]]+theme(strip.text.x=element_text(size=6.8,angle=90),strip.text.y=element_text(size=6.8))
p[[2]][[2]]=p[[2]][[2]]+theme(axis.text.x=element_text(size=6.2))+coord_cartesian(ylim=c(1,15.5))
sy=p[[2]][[2]]$scales$get_scales('y')$clone();sy$limits=c(1,15.5);sy$breaks=c(1,7,13)
p[[2]][[2]]=p[[2]][[2]]+sy
for(j in seq_along(p[[2]][[2]]$layers))if(inherits(p[[2]][[2]]$layers[[j]]$geom,'GeomText')){
 p[[2]][[2]]$layers[[j]]$aes_params$size=6.2/ggplot2::.pt
 p[[2]][[2]]$layers[[j]]$position=position_nudge(y=2)
}
for(i in 1:2)for(j in seq_len(length(p[[i]])))check(old[[i]][[j]],p[[i]][[j]])
save(p,'Figure2',540,565)

for(n in c('b','c')){
 p=readRDS(paste0('results/revision_all_figure_readability_20261109/Figure4',n,'.rds'));old=p
 for(k in c('axis.text.x','axis.text.y','strip.text','legend.text','legend.title')){
  e=ggplot2::calc_element(k,p$theme)
  if(!inherits(e,'element_blank')){e$size=if(k=='legend.title')6.5 else if(n=='b'&&k=='axis.text.y')6.1 else 6.2;p=p+do.call(theme,setNames(list(e),k))}
 }
 p=p+theme(plot.tag.position=c(.015,.98),plot.tag=element_text(size=9,face='bold',vjust=1),plot.margin=margin(4,3,3,3))
 if(n=='b'){
  for(s in p$scales$scales){cat('b scale',s$aesthetics,as.character(s$name),'\n')}
  p=p+labs(fill='Mean Jaccard\nindex',colour='JI > 0.5')+theme(legend.box.spacing=unit(1,'pt'),legend.margin=margin(0,1,0,1),plot.margin=margin(4,3,1,3))
 }
 if(n=='c'){
  e=ggplot2::calc_element('legend.text',p$theme);e$margin=margin(l=2,r=4)
  p=p+theme(legend.text=e,legend.key.width=unit(8,'pt'),legend.spacing.x=unit(4,'pt'))
 }
 check(old,p);save(p,paste0('Figure4',n),c(b=252,c=126)[n],230)
}
p=readRDS('results/revision_all_figure_readability_20261109/Figure4e.rds');old=p
p[[3]]=p[[3]]+coord_cartesian(xlim=c(1,36))
for(j in 1:3)check(old[[j]],p[[j]])
save(p,'Figure4e',540,136)
p=readRDS('results/revision_all_figure_readability_20261109/Supplemental_Figure_S2.rds');old=p
for(j in 1:2){
 e=ggplot2::calc_element('plot.subtitle',p[[1]][[j]]$theme);e$margin=margin(t=10,b=4)
 p[[1]][[j]]=p[[1]][[j]]+theme(plot.subtitle=e)
 check(old[[1]][[j]],p[[1]][[j]])
}
save(p,'Supplemental_Fig_S2',540,560)
writeLines('Data, mappings, all layer data and layer mappings unchanged. Changes limited to annotation position, coordinate display padding, typography and legend spacing.',file.path(out,'plot_verification.txt'))
