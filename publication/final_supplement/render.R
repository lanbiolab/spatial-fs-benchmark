source('results/revision_all_typography_20261109/inputs/fontify.R')
manifest=jsonlite::fromJSON(file.path(out,'inputs/manifest.json'))
for(i in seq_len(nrow(manifest))){
 n=manifest$name[i]
 if(grepl('S10|_S2.pdf',n))next
 p=readRDS(file.path('results/revision_all_figure_readability_20261109',sub('.pdf$','.rds',sub('Supplemental_Fig_','Supplemental_Figure_',n))))
 set.seed(20261109)
 ggsave(file.path(out,n),fontify(p),width=manifest$width[i]/72,height=manifest$height[i]/72,units='in',device=grDevices::cairo_pdf,family='Arial',bg='white')
}
capture=function(path,keep){
 e=new.env(parent=globalenv())
 e$ggsave=function(filename,plot,...){
  n=sub('Supplemental_Figure_','Supplemental_Fig_',basename(filename))
  if(!n %in% keep)return(invisible(NULL))
  args=list(...);args$filename=file.path(out,n);args$plot=fontify(plot)
  args$device=grDevices::cairo_pdf;args$family='Arial';args$bg='white'
  do.call(ggplot2::ggsave,args)
 }
 e$writeLines=function(...)invisible(NULL)
 eval(parse(text=readLines(path)),e)
}
capture(file.path(out,'inputs/render_recheck.R'),'Supplemental_Fig_S2.pdf')
capture(file.path(out,'inputs/render_s10.R'),paste0('Supplemental_Fig_S10',letters[1:6],'.pdf'))
writeLines(audit,file.path(out,'font_data_verification.txt'))
