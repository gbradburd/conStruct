source("~/Dropbox/conStruct/writeup/figs/fig_funcs.R")

setwd("~/Dropbox/conStruct/data/poplars")
n.reps <- 10
K <- 7
for(n in 1:n.reps){
	load(sprintf("poplar_rep%s_test.lnl.Robj",n,n))
	assign(paste0("tl",n),test.lnl)
}
x.vals <- lapply(1:n.reps,function(n){get(sprintf("tl%s",n))})
x.vals.std <- lapply(x.vals,conStruct:::standardize.xvals)
xval.CIs <- conStruct:::get.xval.CIs(x.vals.std,K)
pdf(file="~/Dropbox/conStruct/writeup/figs/populus/populus_std_xval.pdf",width=10,height=5,pointsize=14)
#	quartz(width=10,height=5)
	par(mfrow=c(1,2),mar=c(4,5,4,2))
	plot.xval.CIs(xval.CIs,K,jitter=0.1,xlim=c(0.75,7.25))
		legend(x="bottomright",pch=19,col=c("blue","green"),legend=c("spatial","nonspatial"))
	mtext("Predictive accuracy",side=2,padj=-5,font=2)
	plot.xval.CIs(xval.CIs,K,ylim=c(-40,0),jitter=0.1,xlim=c(1.75,7.25))
		legend(x="bottomright",pch=c(19,NA),lty=c(NA,1),lwd=c(NA,2),col=c(1,adjustcolor(1,0.8)),legend=c("mean","95% CI"))
	mtext("number of clusters",side=1,adj=-1.1,padj=4,font=2)
	mtext("Cross-validation results (Populus)",side=3,adj=5,padj=-2.5,font=2,cex=1.2)
dev.off()


library(maps)
load("~/Dropbox/conStruct/data/poplars/poplar.data.redux.Robj")
col.mat1 <- matrix(ifelse(poplar.data$sp.ID=="Populus trichocarpa","forestgreen","black"),byrow=TRUE,length(poplar.data$sp.ID),length(poplar.data$sp.ID))
col.mat2 <- matrix(ifelse(poplar.data$sp.ID=="Populus trichocarpa","forestgreen","black"),byrow=FALSE,length(poplar.data$sp.ID),length(poplar.data$sp.ID))
pdf(file="~/Dropbox/conStruct/writeup/figs/populus/populus_sampling_map.pdf",width=6,height=5,pointsize=13)
	par(mar=c(4,4,1,1))
	map(xlim = range(poplar.data$coords[,1]) + c(-5,5), ylim = range(poplar.data$coords[,2])+c(-2,2), col="gray")
	points(poplar.data$coords,pch=19,cex=1.2,col=ifelse(poplar.data$sp.ID=="Populus trichocarpa","forestgreen","black"))
	legend(x="bottomleft",pch=19,col=c("forestgreen","black"),legend=c("P. trichocarpa","P. balsamifera"))
	map.axes()
dev.off()


freq.data <- conStruct:::process.freq.data(poplar.data$freqs)
data.block <- conStruct:::make.data.block(K = 1,
										  freq.data = freq.data,
										  coords = poplar.data$coords,
										  spatial = TRUE,
										  geoDist = fields::rdist.earth(poplar.data$coords))

output.list.sp <- vector("list",7)
for(k in 1:7){
	load(sprintf("~/Dropbox/conStruct/data/poplars/runs/poplarsK%s_sp_conStruct.results.Robj",k))
	for(j in 1:k){
		names(conStruct.results[[1]]$MAP$cluster.params[[j]])[4] <- "phi"
	}
	output.list.sp[[k]] <- conStruct.results
}

output.list.nsp <- vector("list",7)
for(k in 1:7){
	load(sprintf("~/Dropbox/conStruct/data/poplars/runs/poplarsK%s_nsp_conStruct.results.Robj",k))
	for(j in 1:k){
		names(conStruct.results[[1]]$MAP$cluster.params[[j]])[4] <- "phi"
	}
	output.list.nsp[[k]] <- conStruct.results
}

csr1.order <- NULL
for(k in 2:7){
	data.block$K <- k
	csr <- output.list.sp[[k]][[1]]
	if(k > 2){
		tmp.csr <- output.list.sp[[k-1]][[1]]
		csr1.order <- match.clusters.x.runs(tmp.csr$MAP$admix.proportions,csr$MAP$admix.proportions,csr1.order)
	}
	if(is.null(csr1.order)){
		csr1.order <- 1:k
	}
	pdf(file=paste0("~/Dropbox/conStruct/writeup/figs/populus/populus_sp",k,".pdf"),width=5,height=5)
		par(mar=c(5,5,1,1))
		map(xlim = range(poplar.data$coords[,1]) + c(-5,5), ylim = range(poplar.data$coords[,2])+c(-2,2), col="gray")
		map.axes()
		make.admix.pie.plot(csr$MAP$admix.proportions[,csr1.order],data.block$coords,cluster.colors=cluster.colors,radii=2.5,add=TRUE)
		box(lwd=2)
	dev.off()
}

pdf(file=paste0("~/Dropbox/conStruct/writeup/figs/populus/populus_sp_clst_covs.pdf"),width=12,height=8,pointsize=14)
	#quartz(width=12,height=8)
	layout(matrix(c(1:6),nrow=2,ncol=3,byrow=TRUE))
	par(mar=c(4,5,3,2),oma=c(3,3,3,1))	
	plot.K.cluster.curves(K=1:6,data.block=data.block,output.list=output.list.sp,col.mat1=col.mat1,col.mat2=col.mat2)
	par(xpd=NA)
		legend(-0.5,0.06,pch=21,
				pt.bg=c(1,"forestgreen","forestgreen"),
				col=c(1,1,"forestgreen"),
				legend=c("balsamifera - balsamifera",
						 "balsamifera - trichocarpa",
						 "trichocarpa - trichocarpa"),cex=0.9,pt.cex=1.5)
	legend(-11,0.057,pch=c(19,NA),lty=c(NA,1),lwd=c(NA,4),legend=c("sample covariance","cluster covariance"))
	mtext(text="geographic distance",side=1,font=2,cex.axis=2,padj=4.5,adj=-3.55)
	mtext(text="allele frequency covariance",side=2,font=2,cex.axis=2,padj=-59.5,adj=20)
dev.off()

pdf(file="~/Dropbox/conStruct/writeup/figs/populus/Fig4_pop_sp_results.pdf",width=12,height=4,pointsize=14)
	csr1.order <- NULL
	par(mfrow=c(1,3),mar=c(5,5,1,1))
		for(k in 2:4){
			csr <- output.list.sp[[k]][[1]]
			if(k > 2){
				tmp.csr <- output.list.sp[[k-1]][[1]]
				csr1.order <- match.clusters.x.runs(tmp.csr$MAP$admix.proportions,csr$MAP$admix.proportions,csr1.order)
			}
			if(is.null(csr1.order)){
				csr1.order <- 1:k
			}
			map(xlim = range(poplar.data$coords[,1]) + c(-5,5), ylim = range(poplar.data$coords[,2])+c(-2,2), col="gray")
			map.axes()
			make.admix.pie.plot(output.list.sp[[k]][[1]]$MAP$admix.proportions[,csr1.order],
								data.block$coords,cluster.colors=cluster.colors,radii=1.7,add=TRUE)
			box(lwd=2)
			mtext(side=1,text=bquote(paste("(",.(letters[k-1]),") ",italic("K")," = ",.(k))),padj=2.7,adj=0.4)
		}
dev.off()

pdf(file="~/Dropbox/conStruct/writeup/figs/populus/populus_nsp_pies.pdf",width=12,height=8,pointsize=14)
	par(mfrow=c(2,3))
	plot.poplar.pies.multipanel(data.block = data.block,
							 K = 7,
							 output.list = output.list.nsp,
							 radii = 1.7,
							 mar = c(5,3,0,1.5))
dev.off()

pdf(file="~/Dropbox/conStruct/writeup/figs/populus/populus_sp_pies.pdf",width=12,height=8,pointsize=14)
	par(mfrow=c(2,3))
	plot.poplar.pies.multipanel(data.block = data.block,
								 K = 7,
								 output.list = output.list.sp,
								 radii = 1.7,
								 mar = c(5,3,0,1.5))
dev.off()

csr1.order <- NULL
for(k in 2:7){
	data.block$K <- k
	data.block$spatial <- FALSE
	csr <- output.list.nsp[[k]][[1]]
	if(k==2){
		tmp.csr <- output.list.sp[[k]][[1]]
		csr1.order <- match.clusters.x.runs(tmp.csr$MAP$admix.proportions,csr$MAP$admix.proportions,csr1.order)
	}
	if(k > 2){
		tmp.csr <- output.list.nsp[[k-1]][[1]]
		csr1.order <- match.clusters.x.runs(tmp.csr$MAP$admix.proportions,csr$MAP$admix.proportions,csr1.order)
	}
	if(is.null(csr1.order)){
		csr1.order <- 1:k
	}
	pdf(file=paste0("~/Dropbox/conStruct/writeup/figs/populus/populus_nsp",k,".pdf"),width=5,height=5)
		par(mar=c(5,5,1,1))
		map(xlim = range(poplar.data$coords[,1]) + c(-5,5), ylim = range(poplar.data$coords[,2])+c(-2,2), col="gray")
		map.axes()
		make.admix.pie.plot(csr$MAP$admix.proportions[,csr1.order],data.block$coords,cluster.colors=cluster.colors,radii=2.5,add=TRUE)
		box(lwd=2)
	dev.off()
}

pdf(file=paste0("~/Dropbox/conStruct/writeup/figs/populus/populus_nsp_clst_covs.pdf"),width=12,height=8,pointsize=14)
	#quartz(width=12,height=8)
	layout(matrix(c(1:6),nrow=2,ncol=3,byrow=TRUE))
	par(mar=c(4,5,3,2),oma=c(3,3,3,1))	
	plot.K.cluster.curves(K=1:6,data.block=data.block,output.list=output.list.nsp,col.mat1=col.mat1,col.mat2=col.mat2,output.list.sp=output.list.sp)
	par(xpd=NA)
		legend(-0.5,0.06,pch=21,
				pt.bg=c(1,"forestgreen","forestgreen"),
				col=c(1,1,"forestgreen"),
				legend=c("balsamifera - balsamifera",
						 "balsamifera - trichocarpa",
						 "trichocarpa - trichocarpa"),cex=0.9,pt.cex=1.5)
	legend(-11,0.057,pch=c(19,NA),lty=c(NA,1),lwd=c(NA,4),legend=c("sample covariance","cluster covariance"))
	mtext(text="geographic distance",side=1,font=2,cex.axis=2,padj=4.5,adj=-3.55)
	mtext(text="allele frequency covariance",side=2,font=2,cex.axis=2,padj=-59.5,adj=20)
dev.off()


K <- 7
laycon.sp <- matrix(0,nrow=7,ncol=7)
colnames(laycon.sp) <- paste(1:7)
for(k in 1:K){
	csr <- output.list.sp[[k]][[1]]
	if(k < 3){
		csr1.order <- NULL
	}
	if(k > 2){
		tmp.csr <- output.list.sp[[k-1]][[1]]
		csr1.order <- match.clusters.x.runs(tmp.csr$MAP$admix.proportions,csr$MAP$admix.proportions,csr1.order)
	}
	if(is.null(csr1.order)){
		csr1.order <- 1:k
	}
	laycon.sp[1:k,k] <- calculate.layer.importance(csr,data.block,csr1.order)
}

laycon.nsp <- matrix(0,nrow=7,ncol=7)
colnames(laycon.nsp) <- paste(1:7)
for(k in 1:K){
	csr <- output.list.nsp[[k]][[1]]
	if(k <= 2){
		csr1.order <- NULL
	}
	if(k==2){
		tmp.csr <- output.list.sp[[k]][[1]]
		csr1.order <- match.clusters.x.runs(tmp.csr$MAP$admix.proportions,csr$MAP$admix.proportions,csr1.order)
	}
	if(k > 2){
		tmp.csr <- output.list.nsp[[k-1]][[1]]
		csr1.order <- match.clusters.x.runs(tmp.csr$MAP$admix.proportions,csr$MAP$admix.proportions,csr1.order)
	}
	if(is.null(csr1.order)){
		csr1.order <- 1:k
	}
	laycon.nsp[1:k,k] <- calculate.layer.importance(csr,data.block,csr1.order)
}

pdf(file="~/Dropbox/conStruct/writeup/figs/populus/populus_laycon_barplots.pdf",width=8,height=4,pointsize=14)
	par(mfrow=c(1,2),mar=c(4,4,3,0.5))
	barplot(laycon.sp,	
			col=cluster.colors,
			xlab="",ylab="layer importance")
	barplot(laycon.nsp,
			col=cluster.colors,
			xlab="",ylab="")
			mtext(side=1,text="number of layers",padj=4.25,adj=-1)
			mtext(side=3,text="Layer contributions (Poplars)",padj=-2.25,adj=7,font=2,cex=1.2)
dev.off()


collapse.ind.Q <- function(Q,pop.vec){
	n.pops <- length(unique(pop.vec))
	admix.props <- matrix(NA,n.pops,ncol(Q))
	for(i in 1:n.pops){
		admix.props[i,] <- colMeans(Q[which(pop.vec==i),,drop=FALSE])
	}
	return(admix.props)
}
load("~/Dropbox/conStruct/data/poplars/fastStructure/poplar.sample.sizes.Robj")
poplar.pop.vec <- unlist(lapply(1:length(poplar.sample.sizes),function(n){rep(n,poplar.sample.sizes[n])}))

pdf(file="~/Dropbox/conStruct/writeup/figs/populus/poplar_fastStr_results.pdf",width=12,height=8,pointsize=14)
	par(mfrow=c(2,3),mar=c(5,3,0,1.5))
	for(k in 2:7){
		w <- collapse.ind.Q(Q = as.matrix(read.table(sprintf("~/Dropbox/conStruct/data/poplars/fastStructure/poplars_K%s.%s.meanQ",k,k),stringsAsFactors=FALSE)),
							pop.vec = poplar.pop.vec)	
		if(k <= 2){
			csr1.order <- NULL
		}
		if(k > 2){
		tmp.w <- collapse.ind.Q(Q = as.matrix(read.table(sprintf("~/Dropbox/conStruct/data/poplars/fastStructure/poplars_K%s.%s.meanQ",k-1,k-1),stringsAsFactors=FALSE)),
								pop.vec = poplar.pop.vec)
			csr1.order <- match.clusters.x.runs(tmp.w,w,csr1.order)
		}
		if(is.null(csr1.order)){
			csr1.order <- 1:k
		}
		par(mar=c(5,5,1,1))
		map(xlim = range(poplar.data$coords[,1]) + c(-5,5), ylim = range(poplar.data$coords[,2])+c(-2,2), col="gray")
		map.axes()
		make.admix.pie.plot(w[,csr1.order],data.block$coords,cluster.colors=cluster.colors,radii=1.7,add=TRUE)
		mtext(side=1,text=bquote(paste("(",.(letters[k-1]),") ",italic("K")," = ",.(k))),padj=2.7,adj=0.4)
	}
dev.off()