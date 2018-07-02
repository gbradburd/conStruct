source("fig_funcs.R")


n.reps <- 10
K <- 7
for(n in 1:n.reps){
	load(sprintf("../../data/poplars/conStruct/xvalidation/poplar_rep%s_test.lnl.Robj",n,n))
	assign(paste0("tl",n),test.lnl)
}
x.vals <- lapply(1:n.reps,function(n){get(sprintf("tl%s",n))})
x.vals.std <- lapply(x.vals,conStruct:::standardize.xvals)
xval.CIs <- conStruct:::get.xval.CIs(x.vals.std,K)
pdf(file="populus/populus_std_xval.pdf",width=10,height=5,pointsize=14)
#	quartz(width=10,height=5)
	par(mfrow=c(1,2),mar=c(4,5,4,2))
	plot.xval.CIs(xval.CIs,K,jitter=0.1,xlim=c(0.75,7.25))
		legend(x="bottomright",pch=19,col=c("blue","green"),legend=c("spatial","nonspatial"))
	mtext("Predictive accuracy",side=2,padj=-5,font=2)
	plot.xval.CIs(xval.CIs,K,ylim=c(-40,0),jitter=0.1,xlim=c(1.75,7.25))
		legend(x="bottomright",pch=c(19,NA),lty=c(NA,1),lwd=c(NA,2),col=c(1,adjustcolor(1,0.8)),legend=c("mean","95% CI"))
	mtext("number of layers",side=1,adj=-0.9,padj=4,font=2)
	mtext("Cross-validation results (Populus)",side=3,adj=5,padj=-2.5,font=2,cex=1.2)
dev.off()


library(maps)
load("../../data/poplars/data/poplar.data.Robj")
col.mat1 <- matrix(ifelse(poplar.data$sp.ID=="Populus trichocarpa","forestgreen","black"),byrow=TRUE,length(poplar.data$sp.ID),length(poplar.data$sp.ID))
col.mat2 <- matrix(ifelse(poplar.data$sp.ID=="Populus trichocarpa","forestgreen","black"),byrow=FALSE,length(poplar.data$sp.ID),length(poplar.data$sp.ID))
pdf(file="populus/populus_sampling_map.pdf",width=6,height=5,pointsize=13)
	par(mar=c(4,4,1,1))
	map(xlim = range(poplar.data$coords[,1]) + c(-5,5), ylim = range(poplar.data$coords[,2])+c(-2,2), col="gray",lforce="e")
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
	load(sprintf("../../data/poplars/conStruct/runs/poplarsK%s_sp_conStruct.results.Robj",k))
	conStruct.results <- cluster.2.layer(conStruct.results)
	for(j in 1:k){
		names(conStruct.results[[1]]$MAP$layer.params[[j]])[4] <- "phi"
	}
	output.list.sp[[k]] <- conStruct.results
}

output.list.nsp <- vector("list",7)
for(k in 1:7){
	load(sprintf("../../data/poplars/conStruct/runs/poplarsK%s_nsp_conStruct.results.Robj",k))
	conStruct.results <- cluster.2.layer(conStruct.results)
	for(j in 1:k){
		names(conStruct.results[[1]]$MAP$layer.params[[j]])[4] <- "phi"
	}
	output.list.nsp[[k]] <- conStruct.results
}

pdf(file=paste0("populus/populus_sp_layer_covs.pdf"),width=12,height=8,pointsize=14)
	#quartz(width=12,height=8)
	layout(matrix(c(1:6),nrow=2,ncol=3,byrow=TRUE))
	par(mar=c(4,5,3,2),oma=c(3,3,3,1))	
	plot.K.layer.curves(K=1:6,data.block=data.block,output.list=output.list.sp,col.mat1=col.mat1,col.mat2=col.mat2)
	par(xpd=NA)
		legend(-0.5,0.06,pch=21,
				pt.bg=c(1,"forestgreen","forestgreen"),
				col=c(1,1,"forestgreen"),
				legend=c("balsamifera - balsamifera",
						 "balsamifera - trichocarpa",
						 "trichocarpa - trichocarpa"),cex=0.9,pt.cex=1.5)
	legend(-11,0.057,pch=c(19,NA),lty=c(NA,1),lwd=c(NA,4),legend=c("sample covariance","layer covariance"))
	mtext(text="geographic distance (mi)",side=1,font=2,cex.axis=2,padj=4.5,adj=-3.55)
	mtext(text="allele frequency covariance",side=2,font=2,cex.axis=2,padj=-59.5,adj=20)
dev.off()

pdf(file=paste0("populus/populus_sp_layer_covs.pdf"),width=12,height=8,pointsize=14)
	#quartz(width=12,height=8)
	layout(matrix(c(1:6),nrow=2,ncol=3,byrow=TRUE))
	par(mar=c(4,5,3,2),oma=c(3,3,3,1))	
	plot.K.layer.curves(K=1:6,data.block=data.block,output.list=output.list.sp,col.mat1=col.mat1,col.mat2=col.mat2)
	par(xpd=NA)
		legend(-0.5,0.06,pch=21,
				pt.bg=c(1,"forestgreen","forestgreen"),
				col=c(1,1,"forestgreen"),
				legend=c("balsamifera - balsamifera",
						 "balsamifera - trichocarpa",
						 "trichocarpa - trichocarpa"),cex=0.9,pt.cex=1.5)
	legend(-11,0.057,pch=c(19,NA),lty=c(NA,1),lwd=c(NA,4),legend=c("sample covariance","layer covariance"))
	mtext(text="geographic distance (mi)",side=1,font=2,cex.axis=2,padj=4.5,adj=-6.2)
	mtext(text="allele frequency covariance",side=2,font=2,cex.axis=2,padj=-59.5,adj=20)
dev.off()

pdf(file="populus/pop_sp_results.pdf",width=12,height=8,pointsize=14)
	#quartz(width=12,height=8,pointsize=14)
	layout(matrix(c(1:6),nrow=2,ncol=3,byrow=TRUE))
	csr1.order <- NULL
	mar <- c(0,3,0,0)
	par(mar=mar,oma=c(4,4,1,1))
		for(k in 2:4){
			csr <- output.list.sp[[k]][[1]]
			if(k > 2){
				tmp.csr <- output.list.sp[[k-1]][[1]]
				csr1.order <- match.layers.x.runs(tmp.csr$MAP$admix.proportions,csr$MAP$admix.proportions,csr1.order)
			}
			if(is.null(csr1.order)){
				csr1.order <- 1:k
			}
			map(xlim = range(poplar.data$coords[,1]) + c(-5,5), ylim = range(poplar.data$coords[,2])+c(-2,2), col="gray",mar=mar,lforce="e")
			map.axes()
			make.admix.pie.plot(output.list.sp[[k]][[1]]$MAP$admix.proportions[,csr1.order],
								data.block$coords,layer.colors=layer.colors,radii=1.7,add=TRUE)
			box(lwd=2)
		mtext(side=1,text=bquote(paste("(",.(letters[k-1]),") ",italic("K")," = ",.(k))),padj=2.7,adj=0.4)
		}
	csr1.order <- NULL
	for(k in 2:4){
		csr <- output.list.sp[[k]][[1]]
		if(k > 2){
			tmp.csr <- output.list.sp[[k-1]][[1]]
			csr1.order <- match.layers.x.runs(tmp.csr$MAP$admix.proportions,csr$MAP$admix.proportions,csr1.order)
		}
		if(is.null(csr1.order)){
			csr1.order <- 1:k
		}
		data.block$K <- k
		plot.layer.curves(data.block = data.block, conStruct.results = output.list.sp[[k]][[1]], layer.cols=layer.colors[order(csr1.order)],add=FALSE,col.mat1=col.mat1,col.mat2=col.mat2)
		box(lwd=2)
		if(k==2){
			mtext(side=2,text="allelic covariance",padj=-3)
			legend(x="topright",pch=c(19,NA),lty=c(NA,1),lwd=c(NA,4),legend=c("sample covariance","layer covariance"))	
		}
		if(k==3){
			mtext(side=1,text="pairwise geographic distance (mi)",padj=3)
		}
		if(k==4){
			legend(x="topright",pch=21,
					pt.bg=c(1,"forestgreen","forestgreen"),
					col=c(1,1,"forestgreen"),
					legend=c("balsamifera - balsamifera",
							 "balsamifera - trichocarpa",
							 "trichocarpa - trichocarpa"),cex=0.9,pt.cex=1.5)
		}
		mtext(side=1,text=bquote(paste("(",.(letters[k+2]),") ",italic("K")," = ",.(k))),padj=4.6,adj=0.4)
	}
dev.off()

pdf(file="populus/populus_nsp_pies.pdf",width=12,height=8,pointsize=14)
	par(mfrow=c(2,3))
	plot.poplar.pies.multipanel(data.block = data.block,
							 K = 7,
							 output.list = output.list.nsp,
							 radii = 1.7,
							 mar = c(5,3,0,1.5),
							 K2.order=c(2,1))
dev.off()

pdf(file="populus/populus_sp_pies.pdf",width=12,height=8,pointsize=14)
	par(mfrow=c(2,3))
	plot.poplar.pies.multipanel(data.block = data.block,
								 K = 7,
								 output.list = output.list.sp,
								 radii = 1.7,
								 mar = c(5,3,0,1.5))
dev.off()

load("../../data/poplars/data/poplar.sample.sizes.Robj")
poplar.pop.vec <- unlist(lapply(1:length(poplar.sample.sizes),function(n){rep(n,poplar.sample.sizes[n])}))

collapse.ind.Q <- function(Q,pop.vec){
	n.pops <- length(unique(pop.vec))
	admix.props <- matrix(NA,n.pops,ncol(Q))
	for(i in 1:n.pops){
		admix.props[i,] <- colMeans(Q[which(pop.vec==i),,drop=FALSE])
	}
	return(admix.props)
}

pdf(file="populus/poplar_admixture_results.pdf",width=12,height=8,pointsize=14)
	par(mfrow=c(2,3),mar=c(5,3,0,1.5))
	for(k in 2:7){
		w <- collapse.ind.Q(Q = as.matrix(read.table(sprintf("../../data/poplars/admixture/poplar.%s.Q",k),stringsAsFactors=FALSE)),
							pop.vec = poplar.pop.vec)
		if(k == 2){
			csr1.order <- c(1,2)
		}
		if(k > 2){
		tmp.w <- collapse.ind.Q(Q = as.matrix(read.table(sprintf("../../data/poplars/admixture/poplar.%s.Q",k-1),stringsAsFactors=FALSE)),
								pop.vec = poplar.pop.vec)
			csr1.order <- match.layers.x.runs(tmp.w,w,csr1.order)
		}
		if(is.null(csr1.order)){
			csr1.order <- 1:k
		}
		par(mar=c(5,5,1,1))
		map(xlim = range(poplar.data$coords[,1]) + c(-5,5), ylim = range(poplar.data$coords[,2])+c(-2,2), col="gray",lforce="e")
		map.axes()
		make.admix.pie.plot(w[,csr1.order],data.block$coords,layer.colors=layer.colors,radii=1.7,add=TRUE)
		mtext(side=1,text=bquote(paste("(",.(letters[k-1]),") ",italic("K")," = ",.(k))),padj=2.7,adj=0.4)
	}
dev.off()

pdf(file="populus/poplars_admixture_CVerror.pdf",width=7,height=7,pointsize=14)
	CV.error <- get.CV.error("../../data/poplars/admixture/exe.admixture.Rout")
	plot(CV.error,pch=19,col="orangered",
			xlab="number of clusters",ylab="Cross-Validation Error",cex=2,
			main="Poplar ADMIXTURE cross-validation results")
	best <- which.min(CV.error)
	points(best,CV.error[best],pch=8,cex=2)
	legend(x="topright",pch=8,legend="model with minimum CV error")
dev.off()