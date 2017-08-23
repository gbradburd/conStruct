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

plot.xval.CIs <- function(xval.CIs,K,k.range=c(1:K),ylim=NULL,cex=1.5,jitter=0,...){
	#recover()
	if(is.null(ylim)){
		ylim <- range(c(unlist(lapply(k.range,function(k){xval.CIs$sp.CIs[[k]]})),
						 unlist(lapply(k.range,function(k){xval.CIs$nsp.CIs[[k]]}))))
	}
	plot(xval.CIs$sp.means,
			ylim=ylim,
			main= "",
			ylab="",
			xlab="",type='n',...)
		lapply(1:K,function(k){
				segments(k,xval.CIs$sp.CIs[[k]][1],
						 k,xval.CIs$sp.CIs[[k]][2],
						 col=adjustcolor(4,0.5),lwd=3)})
		lapply(1:K,function(k){
				segments(k+jitter,xval.CIs$nsp.CIs[[k]][1],
						 k+jitter,xval.CIs$nsp.CIs[[k]][2],
						 col=adjustcolor("green",0.5),lwd=3)})
		points(xval.CIs$sp.means,pch=19,col=4,cex=cex)			
		points(1:K + jitter, xval.CIs$nsp.means,col="green",pch=19,cex=cex)
	return(invisible("plotted"))
}

make.admix.pie.plot.tmp <- function(coords, K, N, conStruct.results, cluster.colors, stat, radii = 2.7, add = FALSE, title = NULL, x.lim = NULL, y.lim = NULL) 
{
    if (is.null(coords)) {
        message("\nuser has not specified sampling coordinates in the data block\n")
    }
    else {
        require(caroline)
        cluster.names <- paste0("cluster_", 1:K)
        sample.names <- paste0("sample_", 1:N)
        color.tab <- nv(c(cluster.colors[1:K]), cluster.names)
        if (stat == "MAP") {
            admix.props <- conStruct.results$MAP$admix.proportions
        }
        else if (stat == "mean") {
            admix.props <- apply(conStruct.results$posterior$admix.proportions, 
                c(2, 3), mean)
        }
        else if (stat == "median") {
            admix.props <- apply(conStruct.results$posterior$admix.proportions, 
                c(2, 3), median)
        }
        pie.list <- lapply(1:N, function(i) {
            nv(admix.props[i, ], cluster.names)
        })
        names(pie.list) <- sample.names
        if (add) {
            par(new = TRUE)
        }
        else {
            par(mar = c(2, 2, 2, 2))
        }
        if (is.null(title)) {
            title <- "Admixture proportion map"
        }
        if (is.null(x.lim)) {
            x.lim <- c(min(coords[, 1]) - 1, max(coords[, 
                1]) + 1)
        }
        if (is.null(y.lim)) {
            y.lim <- c(min(coords[, 2]) - 1, max(coords[, 
                2]) + 1)
        }
        pies(pie.list, x0 = coords[, 1], y0 = coords[, 
            2], color.table = color.tab, border = "black", radii = radii, 
            xlab = "", ylab = "", main = title, lty = 1, density = NULL, 
            xlim = x.lim, ylim = y.lim)
        box(lwd = 2)
    }
    return(invisible(0))
}

assess.significance <- function(mean1,mean2,se1,se2,n){
	sd1 <- se1 * sqrt(n)
	sd2 <- se2 * sqrt(n)
	se <- sqrt( (sd1^2/n + sd2^2/n) )
	df <- ( (sd1^2/n + sd2^2/n)^2 ) / ( (sd1^2/n)^2/(n-1) + (sd2^2/n)^2/(n-1))
	t <- (mean1 - mean2)/se
	p.value <- pt(t,df=df,lower.tail=TRUE)
	return(p.value)
}

choose.best.K <- function(xval.CIs,K,mod,n.reps){
	means <- paste0(mod,".means")
	ses <- paste0(mod,".std.errs")
	overlaps <- lapply(1:(K-1),function(k){
					unlist(
						lapply((k+1):K,function(j){
							assess.significance(xval.CIs[[means]][k],
									xval.CIs[[means]][j],
									xval.CIs[[ses]][k],
									xval.CIs[[ses]][j],
									n.reps)
						})
					)
				})
	best <- FALSE
	k <- 1
	while(!best){
		overlaps
	}
}

choose.best.mod <- function(xval.CIs,K){
	sp.best.K <- choose.best.K(xval.CIs,K,"sp")
}