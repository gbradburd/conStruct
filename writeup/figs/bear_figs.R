source("~/Dropbox/conStruct/writeup/figs/fig_funcs.R")

n.reps <- 10
K <- 7
for(n in 1:n.reps){
	load(sprintf("~/Dropbox/conStruct/data/bears/bear_rep%s_test.lnl.Robj",n,n))
	assign(paste0("tl",n),test.lnl)
}
x.vals <- lapply(1:n.reps,function(n){get(sprintf("tl%s",n))})
x.vals.std <- lapply(x.vals,conStruct:::standardize.xvals)
xval.CIs <- conStruct:::get.xval.CIs(x.vals.std,K)
pdf(file="~/Dropbox/conStruct/writeup/figs/bears/bear_std_xval.pdf",width=10,height=5,pointsize=14)
#	quartz(width=10,height=5)
	par(mfrow=c(1,2),mar=c(4,5,4,2))
	plot.xval.CIs(xval.CIs,K,jitter=0.1,xlim=c(0.75,7.25))
		legend(x="bottomright",pch=19,col=c("blue","green"),legend=c("spatial","nonspatial"))
	mtext("Predictive accuracy",side=2,padj=-5,font=2)
	plot.xval.CIs(xval.CIs,K,ylim=c(-200,0),jitter=0.1,xlim=c(3.75,7.25),xaxt='n')
		axis(1,at=c(4:7))
		legend(x="bottomright",pch=c(19,NA),lty=c(NA,1),lwd=c(NA,2),col=c(1,adjustcolor(1,0.8)),legend=c("mean","95% CI"))
	mtext("number of clusters",side=1,adj=-1.1,padj=4,font=2)
	mtext("Cross-validation results (Bears)",side=3,adj=10.5,padj=-2.5,font=2,cex=1.2)
dev.off()

library(maps)
load("~/Dropbox/conStruct/data/bears/bear.dataset.redux.Robj")
pdf(file="~/Dropbox/conStruct/writeup/figs/bears/bear_sampling_map.pdf",width=6,height=4,pointsize=13)
	#quartz(width=6,height=4,pointsize=13)
	par(mar=c(1,1.5,0,0),oma=c(0,0,0,0))
	map(xlim = range(bear.dataset$sample.coords[,1]) + c(-5,5), ylim = range(bear.dataset$sample.coords[,2])+c(-2,2), col="gray")
	points(bear.dataset$sample.coords,pch=19,cex=1.2)
	map.axes()
dev.off()

freq.data <- conStruct:::process.freq.data(bear.dataset$sample.freqs)
data.block <- conStruct:::make.data.block(K = 1,
										  freq.data = freq.data,
										  coords = bear.dataset$sample.coords,
										  spatial = TRUE,
										  geoDist = fields::rdist.earth(bear.dataset$sample.coords))
output.list.sp <- vector("list",7)
for(k in 1:7){
	load(sprintf("~/Dropbox/conStruct/data/bears/runs/bearsK%s_sp_conStruct.results.Robj",k))
	for(j in 1:k){
		names(conStruct.results[[1]]$MAP$cluster.params[[j]])[4] <- "phi"
	}
	output.list.sp[[k]] <- conStruct.results
}

output.list.nsp <- vector("list",7)
for(k in 1:7){
	load(sprintf("~/Dropbox/conStruct/data/bears/runs/bearsK%s_nsp_conStruct.results.Robj",k))
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
	pdf(file=paste0("~/Dropbox/conStruct/writeup/figs/bears/bear_sp",k,".pdf"),width=7,height=7)
		make.bear.redux.result.plot(admix.proportions = csr$MAP$admix.proportions,
									coords = bear.dataset$sample.coords,
									lump.dist = 200,
									cluster.colors = cluster.colors[order(csr1.order)],
									cluster.order=csr1.order,
									layout = matrix(c(rep(1,10),rep(2,15)),nrow=5,ncol=5,byrow=TRUE))
	dev.off()
}

csr1.order <- NULL
for(k in 2:7){
	data.block$K <- k
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
	pdf(file=paste0("~/Dropbox/conStruct/writeup/figs/bears/bear_nsp",k,".pdf"),width=7,height=7)
		make.bear.redux.result.plot(admix.proportions = csr$MAP$admix.proportions,
									coords = bear.dataset$sample.coords,
									lump.dist = 200,
									cluster.colors = cluster.colors[order(csr1.order)],
									cluster.order=csr1.order,
									layout = matrix(c(rep(1,10),rep(2,15)),nrow=5,ncol=5,byrow=TRUE))
	dev.off()
}

pdf(file=paste0("~/Dropbox/conStruct/writeup/figs/bears/bear_sp_clst_covs.pdf"),width=12,height=8,pointsize=14)
	#quartz(width=12,height=8)
	layout(matrix(c(1:6),nrow=2,ncol=3,byrow=TRUE))
	par(mar=c(4,5,3,2),oma=c(3,3,3,1))	
	plot.K.cluster.curves(K=1:6,data.block=data.block,output.list= output.list.sp,col.mat1=NULL,col.mat2=NULL)
	mtext(text="geographic distance",side=1,font=2,cex.axis=2,padj=4.5,adj=-3.55)
	mtext(text="allele frequency covariance",side=2,font=2,cex.axis=2,padj=-59.5,adj=20)
dev.off()


pdf(file=paste0("~/Dropbox/conStruct/writeup/figs/bears/bear_nsp_clst_covs.pdf"),width=12,height=8,pointsize=14)
	#quartz(width=12,height=8)
	layout(matrix(c(1:6),nrow=2,ncol=3,byrow=TRUE))
	par(mar=c(4,5,3,2),oma=c(3,3,3,1))	
	plot.K.cluster.curves(K=1:6,data.block=data.block,output.list= output.list.nsp,col.mat1=NULL,col.mat2=NULL,output.list.sp=output.list.sp)
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

pdf(file="~/Dropbox/conStruct/writeup/figs/bears/bears_laycon_barplots.pdf",width=8,height=4,pointsize=14)
	par(mfrow=c(1,2),mar=c(4,4,3,0.5))
	barplot(laycon.sp,	
			col=cluster.colors,
			xlab="",ylab="layer importance")
	barplot(laycon.nsp,
			col=cluster.colors,
			xlab="",ylab="")
			mtext(side=1,text="number of layers",padj=4.25,adj=-1)
			mtext(side=3,text="Layer contributions (Bears)",padj=-2.25,adj=18,font=2,cex=1.2)
dev.off()

# pdf(file="~/Dropbox/conStruct/writeup/figs/bears/bears_laycon_sp.pdf",width=7,height=7,pointsize=14)
	# par(mar=c(4.5,4.5,1,1))
	# plot(0,xlim=c(0.7,7.3),ylim=c(0,1),type='n',xlab="number of layers",ylab="layer importance")
		# lapply(1:7,function(k){
			# points(rep(k,k),
				# sort(calculate.layer.contributions(output.list.sp[[k]][[1]]$MAP,data.block),
					# decreasing=TRUE),
				# pch=9,col=c("blue","red","green","yellow","purple","orange","lightblue","darkgreen","lightblue","gray")[1:k],
				# cex=1.5)
		# })
# dev.off()

# pdf(file="~/Dropbox/conStruct/writeup/figs/bears/bears_laycon_nsp.pdf",width=7,height=7,pointsize=14)
	# par(mar=c(4.5,4.5,1,1))
	# plot(0,xlim=c(0.7,7.3),ylim=c(0,1),type='n',xlab="number of layers",ylab="layer importance")
		# lapply(1:7,function(k){
			# points(rep(k,k),
				# sort(calculate.layer.contributions(output.list.nsp[[k]][[1]]$MAP,data.block),
					# decreasing=TRUE),
				# pch=9,col=c("blue","red","green","yellow","purple","orange","lightblue","darkgreen","lightblue","gray")[1:k],
				# cex=1.5)
		# })
# dev.off()


# pdf(file="~/Dropbox/conStruct/writeup/figs/bears/bears_loglaycon_sp.pdf",width=7,height=7,pointsize=14)
	# par(mar=c(4.5,4.5,1,1))
	# plot(0,xlim=c(0.7,7.3),ylim=c(log(0.0008),log(1)),type='n',xlab="number of layers",ylab="log(layer importance)")
		# lapply(1:7,function(k){
			# points(rep(k,k),
				# sort(log(calculate.layer.contributions(output.list.sp[[k]][[1]]$MAP,data.block)),
					# decreasing=TRUE),
				# pch=9,col=c("blue","red","green","yellow","purple","orange","lightblue","darkgreen","lightblue","gray")[1:k],
				# cex=1.5)
		# })
# dev.off()
# pdf(file="~/Dropbox/conStruct/writeup/figs/bears/bears_loglaycon_nsp.pdf",width=7,height=7,pointsize=14)
	# par(mar=c(4.5,4.5,1,1))
	# plot(0,xlim=c(0.7,7.3),ylim=c(log(0.009),log(1)),type='n',xlab="number of layers",ylab="log(layer importance)")
		# lapply(1:7,function(k){
			# points(rep(k,k),
				# sort(log(calculate.layer.contributions(output.list.nsp[[k]][[1]]$MAP,data.block)),
					# decreasing=TRUE),
				# pch=9,col=c("blue","red","green","yellow","purple","orange","lightblue","darkgreen","lightblue","gray")[1:k],
				# cex=1.5)
		# })
# dev.off()
