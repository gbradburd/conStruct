require(conStruct)
cluster.colors <- c("blue","red","goldenrod1","forestgreen","darkorchid1","deepskyblue","darkorange1","seagreen2","yellow1","black")

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

get.admix.CIs <- function(conStruct.results){
	CIs <- apply(conStruct.results$posterior$admix.proportions,
					c(2,3),
					function(x){
						quantile(x,c(0.005,0.995))
					}
			)
	return(CIs)
}

viz.admix.results <- function(sim.admix.props,conStruct.results,clst.order=NULL){
	#layout(matrix(c(1,1,2),1,3))
	k.cols <- c("blue","red","goldenrod1","forestgreen","darkorchid1","deepskyblue","darkorange1","seagreen2","yellow1","black")
	K <- ncol(sim.admix.props)
	if(is.null(clst.order)){
		clst.order <- 1:K
	}
	N <- nrow(sim.admix.props)
	#plot(0,ylim=c(0,1),xlim=c(1,N) + c(-0.5,0.5),type='n',ylab="admixture proportion",xlab="samples")
	admix.prop.CIs <- get.admix.CIs(conStruct.results)
	# for(k in 1:K){
		# segments(1:N,admix.prop.CIs[1,,k],1:N,admix.prop.CIs[2,,k],col=adjustcolor(k.cols[k],0.7),lwd=2.5)
		# points(sim.admix.props[,clst.order[k]],col="gray",bg=k.cols[k],pch=23,cex=1.2)
	# }
	plot(sim.admix.props,conStruct.results$MAP$admix.proportions[,clst.order],type='n',
			xlab="true admixture proportions",
			ylab="estimated admixture proportions",
			main=sprintf("Fitting admixture parameters (K = %s)",K))
		abline(0,1,col=1,lty=2)
		for(k in 1:K){
			segments(sim.admix.props[,clst.order[k]],
						admix.prop.CIs[1,,k],
						sim.admix.props[,clst.order[k]],
						admix.prop.CIs[2,,k],
						col=adjustcolor(k.cols[k],1),lwd=3)
			#points(sim.admix.props[,k],conStruct.results$MAP$admix.proportions[,clst.order[k]],col="gray",bg=k.cols[k],pch=23,cex=1.2)
		}
	legend(x="topleft",lty=1,lwd=2,col=k.cols[1:K],
			legend=paste0("Cluster ",1:K),title="99% credible interval")
	# legend(x="topleft",lty=c(NA,1),lwd=c(NA,2),pch=c(23,NA),col=1,pt.bg="gray",
			# legend=c("true admixture proportion","99% credible interval"))
}

plot.cluster.curves <- function(data.block, conStruct.results, cluster.cols=NULL,sample.cols=NULL,add=FALSE){
	if(is.null(cluster.cols)){
		cluster.colors <- c("blue","red","goldenrod1","forestgreen","darkorchid1","deepskyblue","darkorange1","seagreen2","yellow1","black")
	}
	order.mat <- order(data.block$geoDist)
	if(add==FALSE){
	    y.range <- range(c(
	    			unlist(lapply(1:data.block$K,
	    							function(k){
								        conStruct.results$MAP$cluster.params[[k]]$cluster.cov
					})) + conStruct.results$MAP$gamma, 
					data.block$obsCov))
	    plot(data.block$geoDist[upper.tri(data.block$obsCov, diag = TRUE)], 
	        data.block$obsCov[upper.tri(data.block$obsCov, diag = TRUE)], 
	        xlim = range(data.block$geoDist), ylim = y.range,
	        xlab="",ylab="",pch=19,col=adjustcolor(1,0.7))
    }
    lapply(1:data.block$K, function(k) {
             lines(data.block$geoDist[order.mat],
             		conStruct.results$MAP$gamma + 
             		conStruct.results$MAP$cluster.params[[k]]$cluster.cov[order.mat],
                  col = 1,lwd=4.5,lty=1) ; 
             lines(data.block$geoDist[order.mat],
             		conStruct.results$MAP$gamma + 
             		conStruct.results$MAP$cluster.params[[k]]$cluster.cov[order.mat],
                  col = cluster.cols[k],lwd=4,lty=1)
    })
    return(invisible("cluster covs"))
}

plot.K.cluster.curves <- function(Ks,data.block,output.list,col.mat1=NULL,col.mat2=NULL,output.list.sp=NULL){
cluster.colors <- c("blue","red","goldenrod1","forestgreen","darkorchid1","deepskyblue","darkorange1","seagreen2","yellow1","black")
if(is.null(col.mat1)){
	col.mat1 <- matrix(adjustcolor(1,0.7),data.block$N,data.block$N)
}
if(is.null(col.mat2)){
	col.mat2 <- matrix(adjustcolor(1,0.7),data.block$N,data.block$N)
}
	for(k in Ks){
		data.block$K <- k
		csr <- output.list[[k]][[1]]
		if(k <= 2){
			csr1.order <- NULL
		}
		if(k==2 & !is.null(output.list.sp)){
			tmp.csr <- output.list.sp[[k]][[1]]
			csr1.order <- match.clusters.x.runs(tmp.csr$MAP$admix.proportions,csr$MAP$admix.proportions,csr1.order)
		}
		if(k > 2){
			tmp.csr <- output.list[[k-1]][[1]]
			csr1.order <- match.clusters.x.runs(tmp.csr$MAP$admix.proportions,csr$MAP$admix.proportions,csr1.order)
		}
		if(is.null(csr1.order)){
			csr1.order <- 1:k
		}
				par(mar=c(5,5,1,1))
			    order.mat <- order(data.block$geoDist)
			    y.range <- range(c(
			    			unlist(lapply(1:data.block$K,
			    							function(k){
										        csr$MAP$cluster.params[[k]]$cluster.cov
							})) + csr$MAP$gamma, 
							data.block$obsCov))			
			    plot(data.block$geoDist[upper.tri(data.block$obsCov, diag = TRUE)], 
			        data.block$obsCov[upper.tri(data.block$obsCov, diag = TRUE)], 
			        xlim = range(data.block$geoDist), ylim = y.range,
			        xlab="",ylab="",pch=21,cex=1.3,
			        col=adjustcolor(col.mat1[upper.tri(data.block$obsCov, diag = TRUE)],0.7),
			        bg=adjustcolor(col.mat2[upper.tri(data.block$obsCov, diag = TRUE)],0.7))
				plot.cluster.curves(data.block, csr, cluster.cols=cluster.colors[order(csr1.order)],add=TRUE)
		}
		mtext(text="geographic distance",side=1,font=2,cex.axis=2,padj=4.5,adj=-3.55)
		mtext(text="allele frequency covariance",side=2,font=2,cex.axis=2,padj=-59.5,adj=20)
}

plot.sim.xvals <- function(dir,n.reps,K,simK,y.lim){
	for(n in 1:n.reps){
		load(sprintf("simK%s_rep%s_test.lnl.Robj",simK,n))
		assign(paste0("tl",n),test.lnl)
	}
	x.vals <- lapply(1:n.reps,function(n){get(sprintf("tl%s",n))})
	x.vals.std <- lapply(x.vals,conStruct:::standardize.xvals)
	xval.CIs <- conStruct:::get.xval.CIs(x.vals.std,K)
	plot.xval.CIs(xval.CIs,K,jitter=0.1)
		legend(x="bottomright",pch=19,col=c("blue","green"),legend=c("spatial","nonspatial"))
	mtext("Predictive accuracy",side=2,padj=-5)
	plot.xval.CIs(xval.CIs,K,ylim=y.lim)
		legend(x="bottomleft",pch=c(19,NA),lty=c(NA,1),lwd=c(NA,2),col=c(1,adjustcolor(1,0.8)),legend=c("mean","95% CI"))
	mtext("number of clusters",side=1,adj=-0.95,padj=4)
	mtext(sprintf("Cross-validation results (K=%s)",simK),side=3,adj=70,padj=-2.5,font=2,cex=1.2)
}

plot.sim.pies <- function(data.block,K,output.list,file.name){
	#recover()
	cluster.colors <- c("blue","red","goldenrod1","forestgreen","darkorchid1","deepskyblue","darkorange1","seagreen2","yellow1","black")
	csr1.order <- NULL
	for(k in 2:K){
		data.block$K <- k
		csr <- output.list[[k]][[1]]
		if(k > 2){
			tmp.csr <- output.list[[k-1]][[1]]
			csr1.order <- conStruct:::match.clusters.x.runs(tmp.csr$MAP$admix.proportions,csr$MAP$admix.proportions,csr1.order)
		}
		if(is.null(csr1.order)){
			csr1.order <- 1:k
		}
		pdf(file=paste0(file.name,k,".pdf"),width=5,height=5)
			make.admix.pie.plot(csr$MAP$admix.proportions[,csr1.order],data.block$coords,cluster.colors=cluster.colors,radii=3.5,x.lim=c(2.5,8.5),y.lim=c(2.5,8.5))
		dev.off()
	}
}

# admix.pie.plot <- function (coords, w, cluster.colors, radii = 2.7, add = FALSE,box=TRUE) {
	# K <- ncol(w)
	# N <- nrow(w)
   	# cluster.names <- paste0("cluster_", 1:K)
   	# sample.names <- paste0("sample_", 1:N)
   	# color.tab <- caroline::nv(c(cluster.colors[1:K]),cluster.names)
   	# pie.list <- lapply(1:N, function(i) {
    		   # caroline::nv(w[i, ], cluster.names)
   	# })
   	# names(pie.list) <- sample.names
   	# if (add) {
    	   # par(new = TRUE)
   	# }
   	# else {
    	   # par(mar = c(2, 2, 2, 2))
   	# }
   	# caroline::pies(pie.list, x0 = coords[, 1], 
       # y0 = coords[, 2], color.table = color.tab, 
       # border = "black", radii = radii, xlab = "", ylab = "", 
       # main = "", lty = 1, density = NULL)
	# if(box){
	   	# box(lwd = 2)
	# }
    # return(invisible(0))
# }

# structure.plot <- function(w, mar = c(2, 4, 2, 2), sample.order = NULL, cluster.order = NULL, sample.names = NULL, sort.by = NULL, cluster.colors = NULL) {
	# N <- nrow(w)
	# K <- ncol(w)
    # par(mar = mar)
    # if (is.null(cluster.order)) {
        # cluster.order <- seq(1:K)
    # }
    # if (is.null(sample.order)) {
        # sample.order <- seq(1:N)
    # }
    # if (!is.null(sort.by)) {
        # sample.order <- order(w[,sort.by])
    # }
    # if (is.null(cluster.colors)) {
        # cluster.colors <- c("blue","red","goldenrod1","forestgreen","darkorchid1","deepskyblue","darkorange1","seagreen2","yellow1","black")
    # }
    # use.colors <- cluster.colors[1:K][cluster.order]
    # plot(0, xlim = c(0, N), ylim = c(0, 1), type = "n", 
        # ylab = "", xlab = "", xaxt = "n",bty="n",yaxt="n")
    # plotting.admix.props <- apply(cbind(0, w[, cluster.order]), 1, cumsum)
    # lapply(1:K, function(i) {
        # conStruct:::make.structure.polygon.layer(plotting.admix.props, i, use.colors, sample.order)
    # })
    # if (!is.null(sample.names)) {
        # axis(side = 1, at = seq(1:N) - 0.5, labels = sample.names[sample.order], cex.axis = 0.5, las = 2)
    # }
    # return(invisible("plotted"))
# }

find.unique.sampling.coords <- function(coords,lump.dist){
	sampling.foci <- coords[1,,drop=FALSE]
	focus.membership <- c(1,numeric(nrow(coords)-1))
	for(i in 2:nrow(coords)){
		tmp.dist <- fields::rdist.earth(coords[i,,drop=FALSE],sampling.foci,miles=TRUE)
		if(any(tmp.dist <= lump.dist)){
			focus.membership[i] <- which.min(tmp.dist)
		} else {
			sampling.foci <- rbind(sampling.foci,coords[i,])
			focus.membership[i] <- nrow(sampling.foci)
		}
	}
	return(list("sampling.foci"=sampling.foci,
				"focus.membership"=focus.membership))
}

collapse.rows <- function(matrix,index){
	new.matrix <- matrix(NA,nrow=length(unique(index)),ncol=ncol(matrix))
	for(i in 1:length(unique(index))){
		new.matrix[i,] <- colMeans(matrix[index==i,,drop=FALSE])
	}
	return(new.matrix)
}

divvy.str.by.clump <- function(index,yneg,index.order=NULL){
	#recover()
	N <- length(index)
	n.bins <- length(unique(index))
	if(is.null(index.order)){
		index.order <- 1:n.bins
	}
	bin.size <- sapply(1:n.bins,function(i){length(which(index==i))})
	bin.size <- bin.size[index.order]
	bin.coords <- cbind(c(0,cumsum(bin.size))[1:n.bins],
							cumsum(bin.size))
	bin.middles <- rowMeans(bin.coords)
	segments(x0=bin.coords[,1],y0=0,x1=bin.middles,y1=yneg)
	segments(x0=bin.coords[,2],y0=0,x1=bin.middles,y1=yneg)
	return(cbind(bin.middles,rep(yneg,length(bin.middles))))
}

get.str.plot.order <- function(unique.coords.list){
	n.foci <- nrow(unique.coords.list$sampling.foci)
	n.samples <- length(unique.coords.list$focus.membership)
	ordered.foci <- order(unique.coords.list$sampling.foci[,1])
	str.plot.order <- rep(NA,n.samples)
	ticker <- 1
	for(i in 1:n.foci){
		in.focus <- which(unique.coords.list$focus.membership==ordered.foci[i])
		str.plot.order[in.focus] <- ticker + c(0:(length(in.focus)-1))
		ticker <- max(str.plot.order,na.rm=TRUE) + 1
	}
	return(str.plot.order)
}

make.bear.redux.result.plot <- function(admix.proportions,coords,lump.dist,cluster.colors,cluster.order=NULL,layout=matrix(c(1,2),nrow=2,ncol=1)){
	unique.coords.list <- find.unique.sampling.coords(coords,lump.dist)
	par(xpd=FALSE)
	layout(layout)
	make.structure.plot(admix.proportions = admix.proportions, 
					    mar = c(1,1,1,1),
					    sample.order = order(get.str.plot.order(unique.coords.list)),
					    cluster.order = cluster.order,
					    sample.names = NULL,
					    sort.by = NULL,
					    cluster.colors = cluster.colors)
	bin.middles <- divvy.str.by.clump(unique.coords.list$focus.membership,yneg=-0.03,index.order=order(unique.coords.list$sampling.foci[,1]))
	bin.mid.coords.X <- grconvertX(bin.middles[,1],to="ndc")
	bin.mid.coords.Y <- grconvertY(bin.middles[,2],to="ndc")
	map(xlim = range(coords[,1]) + c(-5,5), ylim = range(coords[,2])+c(-2,2), col="gray",mar=c(1,1,1,1))
	make.admix.pie.plot(coords = unique.coords.list$sampling.foci, 
						admix.proportions = collapse.rows(matrix=admix.proportions,index=unique.coords.list$focus.membership), 
						cluster.colors = cluster.colors,
						radii=3.2,
						add=TRUE)
	box(lwd=2)
	par(xpd=NA)
	segments(x0 = grconvertX(bin.mid.coords.X,from="ndc"),
			 y0 = grconvertY(bin.mid.coords.Y,from="ndc"),
			 x1 = unique.coords.list$sampling.foci[order(unique.coords.list$sampling.foci[,1]),1],
			 y1 = unique.coords.list$sampling.foci[order(unique.coords.list$sampling.foci[,1]),2],
			 lty = 2,
			 lwd=0.5,
			 col = adjustcolor(1,0.6))
}

# match.clusters.x.runs <- function(admix.mat1,admix.mat2,admix.mat1.order=NULL){
	# cluster.colors <- c("blue","red","goldenrod1","forestgreen","darkorchid1","deepskyblue","darkorange1","seagreen2","yellow1","black")
	# K1 <- ncol(run1$MAP$admix.proportions)
		# if(!is.null(run1.order)){
			# run1$MAP$admix.proportions <- run1$MAP$admix.proportions[,run1.order]
		# }
		# K1.cols <- cluster.colors[1:K1]
	# K2 <- ncol(run2$MAP$admix.proportions)
		# K2.cols <- numeric(K2)
	# k.combn <- expand.grid(1:K1,1:K2)
	# clst.sims <- unlist(lapply(1:nrow(k.combn),
						# function(n){
							# measure.frob.similarity(run1$MAP$admix.proportions[,k.combn[n,1],drop=FALSE],
													# run2$MAP$admix.proportions[,k.combn[n,2],drop=FALSE],
													# K=1)
							# }))
	# while(length(which(K2.cols == 0)) > (K2-K1)){
		# tmp.max <- which.max(rank(clst.sims,na.last=FALSE))
		# run2.match <- k.combn[tmp.max,2]
		# run1.match <- k.combn[tmp.max,1]
		# K2.cols[run2.match] <- K1.cols[run1.match]
		# clst.sims[which(k.combn[,1]==run1.match)] <- NA
		# clst.sims[which(k.combn[,2]==run2.match)] <- NA
	# }
	# if(K2 > K1){
		# K2.cols[which(K2.cols==0)] <- cluster.colors[(K1+1):K2]
	# }
	# K2.clst.order <- unique(match(cluster.colors,K2.cols)[which(!is.na(match(cluster.colors,K2.cols)))])
	# clst2.info <- list("cols" = K2.cols,
					   # "clst.order" = K2.clst.order)
	# return(clst2.info)
# }
