#'@export
make.all.the.plots <- function(conStruct.results,n.chains,data.block,prefix,cluster.colors=NULL){
	if(is.null(cluster.colors)){
		cluster.colors <- c("blue","red","green","yellow","purple","orange","lightblue","darkgreen","lightblue","gray")
	}
	lapply(1:n.chains,function(i){
		make.all.chain.plots(conStruct.results[[i]],chain.no=i,data.block,prefix,cluster.colors)
	})
	return(invisible("made chain plots!"))
}

#'@export
make.structure.plot <- function(data.block,conStruct.results,mar=c(2,4,2,2),sample.order=NULL,cluster.order=NULL,sample.names=NULL,sort.by=NULL,cluster.colors=NULL){
	par(mar=mar)
	if(is.null(cluster.order)){
		cluster.order <- seq(1:data.block$K)
	}
	if(is.null(sample.order)){
		sample.order <- seq(1:data.block$N)
	}
	if(!is.null(sort.by)){
		sample.order <- order(conStruct.results$MAP$admix.proportions[,sort.by])
	}
	if(is.null(cluster.colors)){
		cluster.colors <- c("blue","red","green","yellow","purple","orange","lightblue","darkgreen","lightblue","gray")
	}
	if(data.block$K==1){
		conStruct.results$MAP$admix.proportions <- matrix(conStruct.results$MAP$admix.proportions,nrow=data.block$N,ncol=1)
	}
	use.colors <- cluster.colors[1:data.block$K][cluster.order]
	plot(0,xlim=c(0,data.block$N),ylim=c(0,1),type='n',ylab="admixture",xlab="",xaxt='n')
	plotting.admix.props <- apply(cbind(0,conStruct.results$MAP$admix.proportions[,cluster.order]),1,cumsum)
	lapply(1:data.block$K,function(i){
		make.structure.polygon.layer(plotting.admix.props,i,use.colors,sample.order)
	})
	if(!is.null(sample.names)){
		axis(side=1,at=seq(1:data.block$N)-0.5,labels=sample.names[sample.order],cex.axis=0.5,las=2)
	}
	return(invisible("plotted"))
}

#'@export
make.admix.pie.plot <- function(data.block,conStruct.results,cluster.colors,stat,radii=2.7,add=FALSE,title=NULL,x.lim=NULL,y.lim=NULL){
	if(is.null(data.block$coords)){
		message("\nuser has not specified sampling coordinates in the data block\n")
	} else {
		cluster.names <- paste0("cluster_",1:data.block$K)
		sample.names <- paste0("sample_",1:data.block$N)
		color.tab <- caroline::nv(c(cluster.colors[1:data.block$K]),cluster.names)
		if(stat == "MAP"){
			admix.props <- conStruct.results$MAP$admix.proportions
		} else if(stat == "mean"){
			admix.props <- apply(conStruct.results$posterior$admix.proportions,c(2,3),mean)
		} else if(stat == "median"){
			admix.props <- apply(conStruct.results$posterior$admix.proportions,c(2,3),median)		
		}
		pie.list <- lapply(1:data.block$N,function(i){caroline::nv(admix.props[i,],cluster.names)})
		names(pie.list) <- sample.names
		if(add){
			par(new=TRUE)
		} else {
			par(mar=c(2,2,2,2))
		}
		if(is.null(title)){
			title <- "Admixture proportion map"
		}
		if(is.null(x.lim)){
			x.lim <- c(min(data.block$coords[,1]) - 1, max(data.block$coords[,1]) + 1)
		}
		if(is.null(y.lim)){
			y.lim <- c(min(data.block$coords[,2]) - 1, max(data.block$coords[,2]) + 1)
		}
		caroline::pies(pie.list,x0=data.block$coords[,1],y0=data.block$coords[,2],
					color.table=color.tab,border="black",radii=radii,
					xlab="",ylab="",main=title,lty=1,density=NULL,
					xlim = x.lim, ylim = y.lim)
		box(lwd=2)
	}
	return(invisible(0))
}

plot.lpd <- function(conStruct.results){
	plot(conStruct.results$posterior$lpd,
			ylab="posterior probability",
			main="Posterior probability",type='l',
			xlab="MCMC iterations")
	return(invisible(0))
}


plot.nuggets <- function(conStruct.results){
	matplot(conStruct.results$post$nuggets,type='l',
				main="sample nuggets",
				ylab="nugget value",
				xlab="MCMC iterations")
	return(invisible("nuggets"))
}

plot.gamma <- function(conStruct.results){
	plot(conStruct.results$posterior$gamma,
			ylab="gamma",
			xlab="MCMC iterations",
			main="Gamma",type='l')
	return(invisible(0))
}


get.ylim <- function(cluster.params,n.clusters,param){
	y.lim <- range(unlist(
				lapply(
					lapply(1:n.clusters,
						function(i){
							cluster.params[[i]][[param]]
						}),
					function(x){
						range(x)
					})))
	y.lim <- y.lim + c(-0.15*diff(y.lim),0.15*diff(y.lim))
	return(y.lim)
}


plot.cluster.param <- function(cluster.param,clst.col){
	points(cluster.param,type='l',col=clst.col)
	return(invisible(0))
}


plot.cluster.cov.params <- function(data.block,conStruct.results,cluster.colors){
	n.clusters <- data.block$K
	params <- names(conStruct.results$posterior$cluster.params$Cluster_1)[!names(conStruct.results$posterior$cluster.params$Cluster_1)=="cluster.cov"]
	param.ranges <- lapply(params,function(x){get.ylim(conStruct.results$posterior$cluster.params,n.clusters,x)})
	if(length(params) > 0){
		for(i in 1:length(params)){
			plot(0,type='n',main=params[i],
				xlab="MCMC iterations",ylab="parameter value",
				ylim=param.ranges[[i]],xlim=c(1,conStruct.results$posterior$n.iter))
			lapply(1:n.clusters,function(j){plot.cluster.param(conStruct.results$posterior$cluster.params[[j]][[params[i]]],cluster.colors[j])})
			legend(x="topright",col= cluster.colors[1:n.clusters],lty=1,legend=paste0("Cluster_",1:n.clusters))
		}
	}
	return(invisible(0))
}


plot.admix.props <- function(data.block,conStruct.results,cluster.colors){
	n.clusters <- data.block$K
	par(mfrow=c(n.clusters,1),mar=c(3,3,2,2))
		for(i in 1:n.clusters){
			matplot(conStruct.results$posterior$admix.proportions[,,i],type='l',ylim=c(0,1),
					main=paste0("Cluster ",i),ylab="admixture proportion",col=cluster.colors[i])
		}
	return(invisible(0))
}


get.par.cov.CI <- function(data.block,conStruct.results){
	combns <- gtools::combinations(n=data.block$N,r=2,v=1:data.block$N,repeats.allowed=TRUE)
	CIs <- lapply(1:nrow(combns),
				function(i){
					quantile(conStruct.results$posterior$par.cov[,combns[i,1],combns[i,2]],c(0.025,0.975))
				})
	return(CIs)
}

plot.model.fit.CIs <- function(data.block,conStruct.results){
	cov.range <- range(c(data.block$obsCov,
						conStruct.results$posterior$par.cov))
	plot(data.block$geoDist,data.block$obsCov,
    	xlab = "geographic distance", 
        ylab = "covariance",
        main="Cov/geoDist",
        ylim = cov.range, type = "n")
	combns <- gtools::combinations(n=data.block$N,r=2,v=1:data.block$N,repeats.allowed=TRUE)
	CIs <- get.par.cov.CI(data.block,conStruct.results)
	lapply(1:nrow(combns),
			function(i){
				segments(x0 = data.block$geoDist[combns[i,1],combns[i,2]],
						 y0 = CIs[[i]][1],
						 x1 = data.block$geoDist[combns[i,1],combns[i,2]],
						 y1 = CIs[[i]][2],
						 col = adjustcolor(1,0.1),
						 lwd=1.5)
			})
	points(data.block$geoDist,data.block$obsCov,col=2,pch=20,cex=0.8)
	legend(x="topright",legend=c("observed","95% CI"),pch=c(19,NA),lty=c(NA,1),col=c(2,"gray"))
	return(invisible("plotted"))
}


plot.model.fit <- function(data.block,conStruct.results){
	index.mat <- upper.tri(data.block$geoDist, diag = TRUE)
	cov.range <- range(c(data.block$obsCov,conStruct.results$posterior$par.cov))
    plot(data.block$geoDist,data.block$obsCov,
    	xlab = "geographic distance", 
        ylab = "covariance",
        main="Cov/geoDist",
        ylim = cov.range, type = "n")
    lapply(seq(1,conStruct.results$posterior$n.iter,length.out=25), function(i) {
        points(data.block$geoDist[index.mat], conStruct.results$posterior$par.cov[i,,][index.mat],
        	pch = 20, col = adjustcolor(1, 0.1))
    		})
    points(data.block$geoDist[index.mat], data.block$obsCov[index.mat], 
        xlab = "geographic distance", ylab = "covariance", ylim = cov.range, 
        col=2,pch = 19)
	legend(x="topright",legend=c("observed","parametric"),pch=19,col=c(2,1))
	return(invisible("plotted"))
}


plot.cluster.covariances <- function(data.block,conStruct.results,cluster.colors){
	ind.mat <- upper.tri(data.block$geoDist,diag=TRUE)
	order.mat <- order(data.block$geoDist)
	    y.range <- range(c(
	    			unlist(lapply(1:data.block$K,
	    							function(k){
								        conStruct.results$MAP$cluster.params[[k]]$cluster.cov
					})) + conStruct.results$MAP$gamma, 
					data.block$obsCov))
	plot(data.block$geoDist[ind.mat],
		 data.block$obsCov[ind.mat],
			xlim=range(data.block$geoDist),ylim=y.range,
			xlab = "geographic distance",
			ylab = "covariance",
			pch=19,col=adjustcolor(1,0.7))
		lapply(1:data.block$K, function(k) {
             lines(data.block$geoDist[order.mat][ind.mat],
             		conStruct.results$MAP$gamma + 
             		conStruct.results$MAP$cluster.params[[k]]$cluster.cov[order.mat][ind.mat],
                  col = 1,lwd=4.5,lty=1) ; 
             lines(data.block$geoDist[order.mat][ind.mat],
             		conStruct.results$MAP$gamma + 
             		conStruct.results$MAP$cluster.params[[k]]$cluster.cov[order.mat][ind.mat],
                  col = cluster.colors[k],lwd=4,lty=1)
        })
		legend(x="topright",col= cluster.colors[1:data.block$K],lty=1,
				legend=paste0("Cluster_",1:data.block$K),cex=0.7)
	return(invisible("cluster covs"))	
}


structure.polygon <- function(plotting.admix.props,i,j,use.colors){
	polygon(x = c(j-1,j,j,j-1),
			y = c(plotting.admix.props[i,j],
					plotting.admix.props[i,j],
					plotting.admix.props[i+1,j],
					plotting.admix.props[i+1,j]),
			col=use.colors[i])
	return(invisible(j))
}


make.structure.polygon.layer <- function(plotting.admix.props,i,use.colors,sample.order){
		lapply(1:ncol(plotting.admix.props),function(j){
			structure.polygon(plotting.admix.props[,sample.order],i,j,use.colors)
		})
	return(invisible(i))
}

make.all.chain.plots <- function(conStruct.results,chain.no,data.block,prefix,cluster.colors){
	pdf(file=paste0(prefix,"_trace.plots.chain_",chain.no,".pdf"))
		plot.lpd(conStruct.results)
		plot.nuggets(conStruct.results)
		plot.gamma(conStruct.results)
		plot.cluster.cov.params(data.block,conStruct.results,cluster.colors)
		if(data.block$K > 1){
			plot.admix.props(data.block,conStruct.results,cluster.colors)
		}
	dev.off()
	pdf(file=paste0(prefix,"_model.fit.chain_",chain.no,".pdf"))
		plot.model.fit(data.block,conStruct.results)
	dev.off()
	pdf(file=paste0(prefix,"_model.fit.CIs.chain_",chain.no,".pdf"))
		plot.model.fit.CIs(data.block,conStruct.results)
	dev.off()
	if(data.block$spatial | data.block$K > 1){
		pdf(file=paste0(prefix,"_cluster.cov.curves.chain_",chain.no,".pdf"),width=5,height=5)
			plot.cluster.covariances(data.block,conStruct.results,cluster.colors)
		dev.off()
	}
	if(data.block$K > 1){
		pdf(file=paste0(prefix,"_pie.map.chain_",chain.no,".pdf"),width=6,height=6)	
			make.admix.pie.plot(data.block,conStruct.results,cluster.colors,stat="MAP",radii=2.7,add=FALSE,title=NULL,x.lim=NULL,y.lim=NULL)
		dev.off()
		pdf(file=paste0(prefix,"_structure.plot.chain_",chain.no,".pdf"),width=10,height=5)
			make.structure.plot(data.block,conStruct.results,mar=c(2,4,2,2),sample.order=NULL,cluster.order=NULL,sample.names=NULL,sort.by=NULL,cluster.colors=cluster.colors)
		dev.off()
	}
	return(invisible("made chain plots!"))
}
