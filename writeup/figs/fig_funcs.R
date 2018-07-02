require(conStruct)
layer.colors <- c("blue","red","goldenrod1","forestgreen","darkorchid1","deepskyblue","darkorange1","seagreen2","yellow1","black")

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

viz.admix.results <- function(sim.admix.props,conStruct.results,layer.order=NULL){
	#layout(matrix(c(1,1,2),1,3))
	k.cols <- c("blue","red","goldenrod1","forestgreen","darkorchid1","deepskyblue","darkorange1","seagreen2","yellow1","black")
	K <- ncol(sim.admix.props)
	if(is.null(layer.order)){
		layer.order <- 1:K
	}
	N <- nrow(sim.admix.props)
	#plot(0,ylim=c(0,1),xlim=c(1,N) + c(-0.5,0.5),type='n',ylab="admixture proportion",xlab="samples")
	admix.prop.CIs <- get.admix.CIs(conStruct.results)
	# for(k in 1:K){
		# segments(1:N,admix.prop.CIs[1,,k],1:N,admix.prop.CIs[2,,k],col=adjustcolor(k.cols[k],0.7),lwd=2.5)
		# points(sim.admix.props[,layer.order[k]],col="gray",bg=k.cols[k],pch=23,cex=1.2)
	# }
	plot(sim.admix.props,conStruct.results$MAP$admix.proportions[,layer.order],type='n',
			xlab="true admixture proportions",
			ylab="estimated admixture proportions",
			main=bquote(paste("Fitting admixture parameters (true ",italic("K")," = ",.(K),")")))
		abline(0,1,col=1,lty=2)
		for(k in 1:K){
			segments(sim.admix.props[,layer.order[k]],
						admix.prop.CIs[1,,k],
						sim.admix.props[,layer.order[k]],
						admix.prop.CIs[2,,k],
						col=adjustcolor(k.cols[k],1),lwd=3)
			#points(sim.admix.props[,k],conStruct.results$MAP$admix.proportions[,layer.order[k]],col="gray",bg=k.cols[k],pch=23,cex=1.2)
		}
	legend(x="topleft",lty=1,lwd=2,col=k.cols[1:K],
			legend=paste0("Layer ",1:K),title="99% credible interval")
	# legend(x="topleft",lty=c(NA,1),lwd=c(NA,2),pch=c(23,NA),col=1,pt.bg="gray",
			# legend=c("true admixture proportion","99% credible interval"))
}

plot.layer.curves <- function(data.block, conStruct.results, layer.cols=NULL,add=FALSE,col.mat1=NULL,col.mat2=NULL){
	if(is.null(layer.cols)){
		layer.cols <- c("blue","red","goldenrod1","forestgreen","darkorchid1","deepskyblue","darkorange1","seagreen2","yellow1","black")
	}	
	if(is.null(col.mat1)){
		col.mat1 <- matrix(adjustcolor(1,0.7),data.block$N,data.block$N)
	}
	if(is.null(col.mat2)){
		col.mat2 <- matrix(adjustcolor(1,0.7),data.block$N,data.block$N)
	}
	order.mat <- order(data.block$geoDist)
	if(add==FALSE){
	    y.range <- range(c(
	    			unlist(lapply(1:data.block$K,
	    							function(k){
								        conStruct.results$MAP$layer.params[[k]]$layer.cov
					})) + conStruct.results$MAP$gamma, 
					data.block$obsCov)) + c(0,0.02)
	    plot(data.block$geoDist[upper.tri(data.block$obsCov, diag = TRUE)], 
	        data.block$obsCov[upper.tri(data.block$obsCov, diag = TRUE)], 
	        xlim = range(data.block$geoDist), ylim = y.range,
	        xlab="",ylab="",pch=21,cex=1.3,
			col=adjustcolor(col.mat1[upper.tri(data.block$obsCov, diag = TRUE)],0.7),
			bg=adjustcolor(col.mat2[upper.tri(data.block$obsCov, diag = TRUE)],0.7),xaxt='n')
	    axis(side=1,at=c(seq(min(data.block$geoDist),max(data.block$geoDist),length.out=5)),
	    		round(seq(min(data.block$geoDist),max(data.block$geoDist),length.out=5) * sd(fields::rdist.earth(data.block$coords)),0))
    }
    lapply(1:data.block$K, function(k) {
             lines(data.block$geoDist[order.mat],
             		conStruct.results$MAP$gamma + 
             		conStruct.results$MAP$layer.params[[k]]$layer.cov[order.mat],
                  col = 1,lwd=4.5,lty=1) ; 
             lines(data.block$geoDist[order.mat],
             		conStruct.results$MAP$gamma + 
             		conStruct.results$MAP$layer.params[[k]]$layer.cov[order.mat],
                  col = layer.cols[k],lwd=4,lty=1)
    })
    return(invisible("layer covs"))
}

plot.K.layer.curves <- function(Ks,data.block,output.list,col.mat1=NULL,col.mat2=NULL,output.list.sp=NULL,y.range=NULL){
layer.colors <- c("blue","red","goldenrod1","forestgreen","darkorchid1","deepskyblue","darkorange1","seagreen2","yellow1","black")
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
			csr1.order <- match.layers.x.runs(tmp.csr$MAP$admix.proportions,csr$MAP$admix.proportions,csr1.order)
		}
		if(k > 2){
			tmp.csr <- output.list[[k-1]][[1]]
			csr1.order <- match.layers.x.runs(tmp.csr$MAP$admix.proportions,csr$MAP$admix.proportions,csr1.order)
		}
		if(is.null(csr1.order)){
			csr1.order <- 1:k
		}
				par(mar=c(5,5,1,1))
			    order.mat <- order(data.block$geoDist)
			    if(is.null(y.range)){
				    y.range <- range(c(
				    			unlist(lapply(1:data.block$K,
				    							function(k){
											        csr$MAP$layer.params[[k]]$layer.cov
								})) + csr$MAP$gamma, 
								data.block$obsCov))
				}
			    plot(data.block$geoDist[upper.tri(data.block$obsCov, diag = TRUE)], 
			        data.block$obsCov[upper.tri(data.block$obsCov, diag = TRUE)], 
			        xlim = range(data.block$geoDist), ylim = y.range,
			        xlab="",ylab="",pch=21,cex=1.3,
			        col=adjustcolor(col.mat1[upper.tri(data.block$obsCov, diag = TRUE)],0.7),
			        bg=adjustcolor(col.mat2[upper.tri(data.block$obsCov, diag = TRUE)],0.7),xaxt='n')
			    axis(side=1,at=c(seq(min(data.block$geoDist),max(data.block$geoDist),length.out=5)),
			    		round(seq(min(data.block$geoDist),max(data.block$geoDist),length.out=5) * sd(fields::rdist.earth(data.block$coords)),0))
				plot.layer.curves(data.block, csr, layer.cols=layer.colors[order(csr1.order)],add=TRUE)
		}
		#mtext(text="geographic distance (mi)",side=1,font=2,cex.axis=2,padj=4.5,adj=-3.55)
		mtext(text="allele frequency covariance",side=2,font=2,cex.axis=2,padj=-59.5,adj=20)
}

make.col.mats <- function(cols){
	N <- length(cols)
	col.mat1 <- matrix(NA,N,N)
	col.mat2 <- matrix(NA,N,N)
	for(i in 1:N){
		for(j in 1:N){
			col.mat1[i,j] <- cols[i]
			col.mat2[i,j] <- cols[j]
		}
	}
	return(list(col.mat1,col.mat2))
}

plot.sim.xvals <- function(dir,n.reps,K,simK,y.lim){
	#recover()
	for(n in 1:n.reps){
		load(paste0(dir,sprintf("/simK%s_rep%s_test.lnl.Robj",simK,n)))
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
	mtext("number of layers",side=1,adj=-0.85,padj=4)
	mtext(bquote(paste("Cross-validation results (true ",italic("K"),"=",.(simK),")")),side=3,adj=9.25,padj=-2,font=2,cex=1.2)
}

plot.sim.pies <- function(data.block,K,output.list,file.name){
	#recover()
	layer.colors <- c("blue","red","goldenrod1","forestgreen","darkorchid1","deepskyblue","darkorange1","seagreen2","yellow1","black")
	csr1.order <- NULL
	for(k in 2:K){
		data.block$K <- k
		csr <- output.list[[k]][[1]]
		if(k > 2){
			tmp.csr <- output.list[[k-1]][[1]]
			csr1.order <- conStruct:::match.layers.x.runs(tmp.csr$MAP$admix.proportions,csr$MAP$admix.proportions,csr1.order)
		}
		if(is.null(csr1.order)){
			csr1.order <- 1:k
		}
		pdf(file=paste0(file.name,k,".pdf"),width=5,height=5)
			make.admix.pie.plot(csr$MAP$admix.proportions[,csr1.order],data.block$coords,layer.colors=layer.colors,radii=3.5,x.lim=c(2.5,8.5),y.lim=c(2.5,8.5))
		dev.off()
	}
}

plot.sim.pies.multipanel <- function(data.block,K,output.list,radii,mar=NULL,trueK){
	#recover()
	layer.colors <- c("blue","red","goldenrod1","forestgreen","darkorchid1","deepskyblue","darkorange1","seagreen2","yellow1","black")
	csr1.order <- NULL
	for(k in 2:K){
		data.block$K <- k
		csr <- output.list[[k]][[1]]
		if(k > 2){
			tmp.csr <- output.list[[k-1]][[1]]
			csr1.order <- conStruct:::match.layers.x.runs(tmp.csr$MAP$admix.proportions,csr$MAP$admix.proportions,csr1.order)
		}
		if(is.null(csr1.order)){
			csr1.order <- 1:k
		}
		make.admix.pie.plot(csr$MAP$admix.proportions[,csr1.order],data.block$coords,layer.colors=layer.colors,radii=radii,x.lim=c(2.5,8.5),y.lim=c(2.5,8.5),mar=mar)
		if(k==3){
			mtext(side=3,text=bquote(paste("True ",italic("K")," = ",.(trueK))),cex=1.5,padj=-0.7)
		}
			mtext(side=1,text=bquote(paste("(",.(letters[k-1]),") ",italic("K")," = ",.(k))),padj=2.7,adj=0.4)
	}
}

plot.poplar.pies.multipanel <- function(data.block,K,output.list,radii,mar=NULL,K2.order=NULL){
	#recover()
	layer.colors <- c("blue","red","goldenrod1","forestgreen","darkorchid1","deepskyblue","darkorange1","seagreen2","yellow1","black")
	csr1.order <- NULL
	for(k in 2:K){
		data.block$K <- k
		csr <- output.list[[k]][[1]]
		if(k==2 & !is.null(K2.order)){
			csr1.order <- K2.order
		}
		if(k > 2){
			tmp.csr <- output.list[[k-1]][[1]]
			csr1.order <- conStruct:::match.layers.x.runs(tmp.csr$MAP$admix.proportions,csr$MAP$admix.proportions,csr1.order)
		}
		if(is.null(csr1.order)){
			csr1.order <- 1:k
		}
		map(xlim = range(poplar.data$coords[,1]) + c(-5,5), ylim = range(poplar.data$coords[,2])+c(-2,2), col="gray", mar=mar)
		map.axes()
		make.admix.pie.plot(csr$MAP$admix.proportions[,csr1.order],data.block$coords,layer.colors=layer.colors,radii=radii,x.lim=c(2.5,8.5),y.lim=c(2.5,8.5),add=TRUE)
			mtext(side=1,text=bquote(paste("(",.(letters[k-1]),") ",italic("K")," = ",.(k))),padj=2.7,adj=0.4)
	}
}

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

make.bear.redux.result.plot <- function(admix.proportions,coords,lump.dist,layer.colors,layer.order=NULL,layout=matrix(c(1,2),nrow=2,ncol=1)){
	unique.coords.list <- find.unique.sampling.coords(coords,lump.dist)
	par(xpd=FALSE)
	layout(layout)
	make.structure.plot(admix.proportions = admix.proportions, 
					    mar = c(1,1,1,1),
					    sample.order = order(get.str.plot.order(unique.coords.list)),
					    layer.order = layer.order,
					    sample.names = NULL,
					    sort.by = NULL,
					    layer.colors = layer.colors)
	bin.middles <- divvy.str.by.clump(unique.coords.list$focus.membership,yneg=-0.03,index.order=order(unique.coords.list$sampling.foci[,1]))
	bin.mid.coords.X <- grconvertX(bin.middles[,1],to="ndc")
	bin.mid.coords.Y <- grconvertY(bin.middles[,2],to="ndc")
	map(xlim = range(coords[,1]) + c(-5,5), ylim = range(coords[,2])+c(-2,2), col="gray",mar=c(1,1,1,1))
	make.admix.pie.plot(coords = unique.coords.list$sampling.foci, 
						admix.proportions = collapse.rows(matrix=admix.proportions,index=unique.coords.list$focus.membership), 
						layer.colors = layer.colors,
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

make.bear.redux.result.plot.multipanel1 <- function(admix.proportions,coords,lump.dist,layer.colors,layer.order=NULL){
	unique.coords.list <- find.unique.sampling.coords(coords,lump.dist)
	par(xpd=FALSE)
	make.structure.plot(admix.proportions = admix.proportions, 
					    mar = c(1,2.5,1,1),
					    sample.order = order(get.str.plot.order(unique.coords.list)),
					    layer.order = layer.order,
					    sample.names = NULL,
					    sort.by = NULL,
					    layer.colors = layer.colors)
	bin.middles <- divvy.str.by.clump(unique.coords.list$focus.membership,yneg=-0.03,index.order=order(unique.coords.list$sampling.foci[,1]))
	bin.mid.coords.X <- grconvertX(bin.middles[,1],to="ndc")
	bin.mid.coords.Y <- grconvertY(bin.middles[,2],to="ndc")
	map(xlim = range(coords[,1]) + c(-5,5), ylim = range(coords[,2])+c(-2,2), col="gray",mar=c(5,0.5,0.5,0.5),lforce="e")
	make.admix.pie.plot(coords = unique.coords.list$sampling.foci, 
						admix.proportions = collapse.rows(matrix=admix.proportions,index=unique.coords.list$focus.membership), 
						layer.colors = layer.colors,
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

make.bear.redux.result.plot.multipanel2 <- function(output.list,coords,lump.dist,layer.colors,csr1.order=NULL){
	unique.coords.list <- find.unique.sampling.coords(coords,lump.dist)
	for(k in 2:7){
		data.block$K <- k
		csr <- output.list[[k]][[1]]
		if(k > 2){
			tmp.csr <- output.list[[k-1]][[1]]
			csr1.order <- match.layers.x.runs(tmp.csr$MAP$admix.proportions,csr$MAP$admix.proportions,csr1.order)
		}
		if(is.null(csr1.order)){
			csr1.order <- 1:k
		}
		par(xpd=FALSE)
		make.structure.plot(admix.proportions = csr$MAP$admix.proportions, 
						    mar = c(1,2.5,1,1),
						    sample.order = order(get.str.plot.order(unique.coords.list)),
						    layer.order = csr1.order,
						    sample.names = NULL,
						    sort.by = NULL,
						    layer.colors = layer.colors[order(csr1.order)])
		bin.middles <- divvy.str.by.clump(unique.coords.list$focus.membership,yneg=-0.03,index.order=order(unique.coords.list$sampling.foci[,1]))
		bin.mid.coords.X <- grconvertX(bin.middles[,1],to="ndc")
		bin.mid.coords.Y <- grconvertY(bin.middles[,2],to="ndc")
		map(xlim = range(coords[,1]) + c(-5,5), ylim = range(coords[,2])+c(-2,2), col="gray",mar=c(5,0.5,0.5,0.5),lforce="e")
		make.admix.pie.plot(coords = unique.coords.list$sampling.foci, 
							admix.proportions = collapse.rows(matrix=csr$MAP$admix.proportions,index=unique.coords.list$focus.membership), 
							layer.colors = layer.colors[order(csr1.order)],
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
		mtext(side=1,text=bquote(paste("(",.(letters[k-1]),") ",italic("K")," = ",.(k))),padj=-1.5,adj=0.03,cex=1.3)
	}
}

cluster.2.layer <- function(conStruct.results){
	names(conStruct.results[[1]][["posterior"]])[which(names(conStruct.results[[1]][[1]])=="cluster.params")] <- "layer.params"
		names(conStruct.results[[1]]$posterior$layer.params) <- gsub("Cluster","Layer",names(conStruct.results[[1]]$posterior$layer.params))
			for(k in 1:length(conStruct.results[[1]]$posterior$layer.params)){
				names(conStruct.results[[1]]$posterior$layer.params[[k]])[which(names(conStruct.results[[1]]$posterior$layer.params[[k]])=="cluster.cov")] <- "layer.cov"
			}
	names(conStruct.results[[1]][["MAP"]])[which(names(conStruct.results[[1]]$MAP)=="cluster.params")] <- "layer.params"
		names(conStruct.results[[1]]$MAP$layer.params) <- gsub("Cluster","Layer",names(conStruct.results[[1]]$MAP$layer.params))
			for(k in 1:length(conStruct.results[[1]]$MAP$layer.params)){
				names(conStruct.results[[1]]$MAP$layer.params[[k]])[which(names(conStruct.results[[1]]$MAP$layer.params[[k]])=="cluster.cov")] <- "layer.cov"
			}
	return(conStruct.results)
}

get.CV.error <- function(Rout.file){
	log <- scan(Rout.file,what="character",sep="\n")
	CV.error <- as.numeric(
					unlist(
						lapply(
							strsplit(
								log[grepl("CV error",log)],
								": "),
						"[[",2)
					)
				)
	return(CV.error)
}