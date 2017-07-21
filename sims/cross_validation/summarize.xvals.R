#summary of cross-validation results

standardize.xvals <- function(x.val){
	mean.lnls <- lapply(x.val,function(x){
					lapply(x,function(k){mean(unlist(k))})
				 })
	xval.max <- max(unlist(mean.lnls))
	mean.std.lnls <- lapply(mean.lnls,function(s){
						lapply(s,function(k){
							k - xval.max
						})})
	return(mean.std.lnls)
}

get.xval.CIs <- function(x.vals.std,K){
	#recover()
	sp.means <- lapply(
				   		lapply(1:K,
							function(k){
								unlist(lapply(
									lapply(x.vals.std,"[[","sp"),
						"[[",k))}),
					function(x){
						mean(x)})
	sp.std.errs <- lapply(
				   		lapply(1:K,
							function(k){
								unlist(lapply(
									lapply(x.vals.std,"[[","sp"),
						"[[",k))}),
					function(x){
						sd(x)/sqrt(length(x))})
	sp.CIs <- lapply(1:K,function(k){
					sp.means[[k]] + c(-1.96*sp.std.errs[[k]],
									   1.96*sp.std.errs[[k]])})
	nsp.means <- lapply(
				   		lapply(1:K,
							function(k){
								unlist(lapply(
									lapply(x.vals.std,"[[","nsp"),
						"[[",k))}),
					function(x){
						mean(x)})
	nsp.std.errs <- lapply(
				   		lapply(1:K,
							function(k){
								unlist(lapply(
									lapply(x.vals.std,"[[","nsp"),
						"[[",k))}),
					function(x){
						sd(x)/sqrt(length(x))})
	nsp.CIs <- lapply(1:K,function(k){
					nsp.means[[k]] + c(-1.96*nsp.std.errs[[k]],
									    1.96*nsp.std.errs[[k]])})
	return(list("sp.means" = unlist(sp.means),
				"sp.std.errs" = unlist(sp.std.errs),
				"sp.CIs" = sp.CIs,
				"nsp.means" = unlist(nsp.means),
				"nsp.std.errs" = unlist(nsp.std.errs),
				"nsp.CIs" = nsp.CIs))
}

plot.xval.CIs <- function(xval.CIs,K,k.range=c(1:K),ylim=NULL,simK=NULL){
	#recover()
	if(is.null(ylim)){
		ylim <- range(c(unlist(lapply(k.range,function(k){xval.CIs$sp.CIs[[k]]})),
						 unlist(lapply(k.range,function(k){xval.CIs$nsp.CIs[[k]]}))))
	}
	plot(xval.CIs$sp.means,
			ylim=ylim,
			xlim=range(k.range),
			main= paste0("Cross-validation results",simK),
			ylab="mean log p(testing data partition | training params)",
			xlab="number of clusters",type='n')
		lapply(1:K,function(k){
				segments(k,xval.CIs$sp.CIs[[k]][1],
						 k,xval.CIs$sp.CIs[[k]][2],
						 col=adjustcolor(4,0.5),lwd=3)})
		lapply(1:K,function(k){
				segments(k,xval.CIs$nsp.CIs[[k]][1],
						 k,xval.CIs$nsp.CIs[[k]][2],
						 col=adjustcolor("green",0.5),lwd=3)})
		points(xval.CIs$sp.means,pch=19,col=4,cex=1.5,)			
		points(xval.CIs$nsp.means,col="green",pch=19,cex=1.5)
	return(invisible("plotted"))
}