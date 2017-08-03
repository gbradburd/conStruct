source("~/Dropbox/conStruct/sims/cross_validation/summarize.xvals.R")

plot.cluster.curves <- function(data.block, conStruct.results, cluster.cols=NULL,sample.cols=NULL,add=FALSE){
	if(is.null(cluster.cols)){
		cluster.cols <- c("blue", "red", "green", "yellow", "purple","orange", "lightblue", "darkgreen", "lightblue", "gray")
	}
	if(add==FALSE){
		order.mat <- order(data.block$geoDist)
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


setwd("~/Dropbox/gid_runs/populus/all/runs")
n.reps <- 10
K <- 7
for(n in 1:n.reps){
	load(sprintf("poplar_rep%srep%s_test.lnl.Robj",n,n))
	assign(paste0("tl",n),test.lnl)
}
x.vals <- lapply(1:n.reps,function(n){get(sprintf("tl%s",n))})
x.vals.std <- lapply(x.vals,standardize.xvals)
xval.CIs <- get.xval.CIs(x.vals.std,K)
pdf(file="~/Dropbox/conStruct/writeup/figs/populus/populus_std_xval.pdf",width=10,height=5,pointsize=14)
#	quartz(width=10,height=5)
	par(mfrow=c(1,2),mar=c(4,5,4,2))
	plot.xval.CIs(xval.CIs,K,jitter=0.1,xlim=c(0.75,7.25))
		legend(x="bottomright",pch=19,col=c("blue","green"),legend=c("spatial","nonspatial"))
	mtext("Predictive accuracy",side=2,padj=-5,font=2)
	plot.xval.CIs(xval.CIs,K,ylim=c(-100,0),jitter=0.1,xlim=c(1.75,7.25))
		legend(x="bottomright",pch=c(19,NA),lty=c(NA,1),lwd=c(NA,2),col=c(1,adjustcolor(1,0.8)),legend=c("mean","95% CI"))
	mtext("number of clusters",side=1,adj=-1.1,padj=4,font=2)
	mtext("Cross-validation results (Populus)",side=3,adj=5,padj=-2.5,font=2,cex=1.2)
dev.off()


library(conStruct)
library(maps)
my.rep <- 1
load("~/Dropbox/gid_runs/populus/poplar.spStr.dataset.Robj")
col.mat1 <- matrix(ifelse(poplar.data$sp.ID=="Populus trichocarpa","forestgreen","black"),byrow=TRUE,length(poplar.data$sp.ID),length(poplar.data$sp.ID))
col.mat2 <- matrix(ifelse(poplar.data$sp.ID=="Populus trichocarpa","forestgreen","black"),byrow=FALSE,length(poplar.data$sp.ID),length(poplar.data$sp.ID))
pdf(file="~/Dropbox/conStruct/writeup/figs/populus/populus_sampling_map.pdf",width=6,height=5,pointsize=13)
	par(mar=c(4,4,1,1))
	map(xlim = range(poplar.data$coords[,1]) + c(-5,5), ylim = range(poplar.data$coords[,2])+c(-2,2), col="gray")
	points(poplar.data$coords,pch=19,cex=1.2,col=ifelse(poplar.data$sp.ID=="Populus trichocarpa","forestgreen","black"))
	legend(x="bottomleft",pch=19,col=c("forestgreen","black"),legend=c("P. trichocarpa","P. balsamifera"))
	map.axes()
dev.off()



for(k in 2:7){
	load(sprintf("~/Dropbox/gid_runs/populus/all/runs/poplar_rep%s_sp_rep%sK%s_freq.data.Robj",my.rep,my.rep,k))
	data.block <- make.data.block(K = k,
								  freq.data = freq.data,
								  coords = poplar.data$coords,
								  spatial = TRUE,
								  geoDist = fields::rdist.earth(poplar.data$coords))
	load(sprintf("~/Dropbox/gid_runs/populus/all/runs/poplar_rep%s_sp_rep%sK%s_conStruct.results.Robj",my.rep,my.rep,k))
	csr <- conStruct.results[[1]]
	clst.match <- NULL
	if(k < 4){
		csr1.order <- NULL
	}
	clst.match$cols <- c(4,2)
	if(k > 2){
	load(sprintf("~/Dropbox/gid_runs/populus/all/runs/poplar_rep%s_sp_rep%sK%s_conStruct.results.Robj",my.rep,my.rep,k-1))
		tmp.csr <- conStruct.results[[1]]
		clst.match <- match.clusters.x.runs(tmp.csr,csr,csr1.order)
		csr1.order <- clst.match$clst.order
	}
	# pdf(file=paste0("~/Dropbox/conStruct/writeup/figs/sims/simK1_nsp_str_K",k,".pdf"),width=10,height=5)
		# make.structure.plot(data.block,csr,cluster.order=clst.match$clst.order,cluster.colors=clst.match$cols)
	# dev.off()
	pdf(file=paste0("~/Dropbox/conStruct/writeup/figs/populus/populus_sp",k,".pdf"),width=5,height=5)
		par(mar=c(5,5,1,1))
		map(xlim = range(poplar.data$coords[,1]) + c(-5,5), ylim = range(poplar.data$coords[,2])+c(-2,2), col="gray")
		map.axes()
		make.admix.pie.plot(data.block,csr,cluster.colors=clst.match$cols,stat="MAP",title="",radii=2.5,add=TRUE)
	dev.off()
}

pdf(file=paste0("~/Dropbox/conStruct/writeup/figs/populus/populus_sp_clst_covs.pdf"),width=12,height=8,pointsize=14)
	#quartz(width=12,height=8)
	layout(matrix(c(1:6),nrow=2,ncol=3,byrow=TRUE))
	par(mar=c(4,5,3,2),oma=c(3,3,3,1))	
	for(k in 2:7){
		load(sprintf("~/Dropbox/gid_runs/populus/all/runs/poplar_rep%s_sp_rep%sK%s_freq.data.Robj",my.rep,my.rep,k))
		data.block <- make.data.block(K = k,
									  freq.data = freq.data,
									  coords = poplar.data$coords,
									  spatial = TRUE,
									  geoDist = fields::rdist.earth(poplar.data$coords))
		load(sprintf("~/Dropbox/gid_runs/populus/all/runs/poplar_rep%s_sp_rep%sK%s_conStruct.results.Robj",my.rep,my.rep,k))
		csr <- conStruct.results[[1]]
		clst.match <- NULL
		if(k < 4){
			csr1.order <- NULL
		}
		clst.match$cols <- c(4,2)
		if(k > 2){
		load(sprintf("~/Dropbox/gid_runs/populus/all/runs/poplar_rep%s_sp_rep%sK%s_conStruct.results.Robj",my.rep,my.rep,k-1))
			tmp.csr <- conStruct.results[[1]]
			clst.match <- match.clusters.x.runs(tmp.csr,csr,csr1.order)
			csr1.order <- clst.match$clst.order
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
			plot.cluster.curves(data.block, csr, cluster.cols=clst.match$cols,add=TRUE)
			if(k==2){
				legend(x="topright",pch=c(19,NA),lty=c(NA,1),lwd=c(NA,4),legend=c("sample covariance","cluster covariance"))
			}
			if(k==7){
				legend(x="topright",pch=21,
						pt.bg=c(1,"forestgreen","forestgreen"),
						col=c(1,1,"forestgreen"),
						legend=c("balsamifera - balsamifera",
								 "balsamifera - trichocarpa",
								 "trichocarpa - trichocarpa"),cex=0.9,pt.cex=1.5)
			}
	}
	mtext(text="geographic distance",side=1,font=2,cex.axis=2,padj=4.5,adj=-3.55)
	mtext(text="allele frequency covariance",side=2,font=2,cex.axis=2,padj=-59.5,adj=20)
dev.off()



for(k in 2:7){
	load(sprintf("~/Dropbox/gid_runs/populus/all/runs/poplar_rep%s_nsp_rep%sK%s_freq.data.Robj",my.rep,my.rep,k))
	data.block <- make.data.block(K = k,
								  freq.data = freq.data,
								  coords = poplar.data$coords,
								  spatial = FALSE,
								  geoDist = fields::rdist.earth(poplar.data$coords))
	load(sprintf("~/Dropbox/gid_runs/populus/all/runs/poplar_rep%s_nsp_rep%sK%s_conStruct.results.Robj",my.rep,my.rep,k))
	csr <- conStruct.results[[1]]
	clst.match <- NULL
	if(k < 4){
		csr1.order <- NULL
	}
	clst.match$cols <- c(4,2)
	if(k > 2){
	load(sprintf("~/Dropbox/gid_runs/populus/all/runs/poplar_rep%s_nsp_rep%sK%s_conStruct.results.Robj",my.rep,my.rep,k-1))
		tmp.csr <- conStruct.results[[1]]
		clst.match <- match.clusters.x.runs(tmp.csr,csr,csr1.order)
		csr1.order <- clst.match$clst.order
	}
	# pdf(file=paste0("~/Dropbox/conStruct/writeup/figs/sims/simK1_nsp_str_K",k,".pdf"),width=10,height=5)
		# make.structure.plot(data.block,csr,cluster.order=clst.match$clst.order,cluster.colors=clst.match$cols)
	# dev.off()
	pdf(file=paste0("~/Dropbox/conStruct/writeup/figs/populus/populus_nsp",k,".pdf"),width=5,height=5)
		map(xlim = range(poplar.data$coords[,1]) + c(-5,5), ylim = range(poplar.data$coords[,2])+c(-2,2), col="gray")
		map.axes()
		make.admix.pie.plot(data.block,csr,cluster.colors=clst.match$cols,stat="MAP",title="",radii=2.5,add=TRUE)
	dev.off()
}

pdf(file=paste0("~/Dropbox/conStruct/writeup/figs/populus/populus_nsp_clst_covs.pdf"),width=12,height=8,pointsize=14)
	#quartz(width=12,height=8)
	layout(matrix(c(1:6),nrow=2,ncol=3,byrow=TRUE))
	par(mar=c(4,5,3,2),oma=c(3,3,3,1))	
	for(k in 2:7){
		load(sprintf("~/Dropbox/gid_runs/populus/all/runs/poplar_rep%s_nsp_rep%sK%s_freq.data.Robj",my.rep,my.rep,k))
		data.block <- make.data.block(K = k,
									  freq.data = freq.data,
									  coords = poplar.data$coords,
									  spatial = FALSE,
									  geoDist = fields::rdist.earth(poplar.data$coords))
		load(sprintf("~/Dropbox/gid_runs/populus/all/runs/poplar_rep%s_nsp_rep%sK%s_conStruct.results.Robj",my.rep,my.rep,k))
		csr <- conStruct.results[[1]]
		clst.match <- NULL
		if(k < 4){
			csr1.order <- NULL
		}
		clst.match$cols <- c(4,2)
		if(k > 2){
		load(sprintf("~/Dropbox/gid_runs/populus/all/runs/poplar_rep%s_nsp_rep%sK%s_conStruct.results.Robj",my.rep,my.rep,k-1))
			tmp.csr <- conStruct.results[[1]]
			clst.match <- match.clusters.x.runs(tmp.csr,csr,csr1.order)
			csr1.order <- clst.match$clst.order
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
			plot.cluster.curves(data.block, csr, cluster.cols=clst.match$cols,add=TRUE)
			if(k==2){
				legend(x="topright",pch=c(19,NA),lty=c(NA,1),lwd=c(NA,4),legend=c("sample covariance","cluster covariance"))
			}
			if(k==5){
				legend(x="topright",pch=21,
						pt.bg=c(1,"forestgreen","forestgreen"),
						col=c(1,1,"forestgreen"),
						legend=c("balsamifera - balsamifera",
								 "balsamifera - trichocarpa",
								 "trichocarpa - trichocarpa"),cex=0.7,pt.cex=1.3)
			}
	}
	mtext(text="geographic distance",side=1,font=2,cex.axis=2,padj=4.5,adj=-3.55)
	mtext(text="allele frequency covariance",side=2,font=2,cex.axis=2,padj=-59.5,adj=20)
dev.off()