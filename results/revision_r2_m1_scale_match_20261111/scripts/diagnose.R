o <- 'results/revision_r2_m1_scale_match_20261111'
d <- read.delim(file.path(o,'analysis/pooled_coordinates.tsv'))
f <- getFromNamespace('transloc_func','SPARK')
rows <- list()
for(sid in unique(d$sample_id)) {
 xy <- as.matrix(d[d$sample_id==sid,c('x','y')])
 for(kernel in c('gaussian','cosine')) for(k in 1:5) {
  z <- apply(xy,2,f,lker=k,transfunc=kernel)
  z <- apply(z,2,scale,scale=FALSE);g<-crossprod(z,z)
  err<-tryCatch({solve(g);''},error=function(e)conditionMessage(e))
  rows[[length(rows)+1L]]<-data.frame(sample_id=sid,kernel=kernel,k=k,x_variance=var(z[,1]),y_variance=var(z[,2]),determinant=det(g),rcond=rcond(g),solve_error=err)
 }
}
write.table(do.call(rbind,rows),file.path(o,'analysis/sparkx_kernel_diagnostics.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
print(subset(do.call(rbind,rows),nzchar(solve_error)))
