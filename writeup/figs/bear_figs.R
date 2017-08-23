source("~/Dropbox/conStruct/writeup/figs/fig_funcs.R")

n.reps <- 6
K <- 7
for(n in 1:n.reps){
	load(sprintf("bear_rep%s_test.lnl.Robj",n,n))
	assign(paste0("tl",n),test.lnl)
}
x.vals <- lapply(1:n.reps,function(n){get(sprintf("tl%s",n))})
x.vals.std <- lapply(x.vals,standardize.xvals)
xval.CIs <- get.xval.CIs(x.vals.std,K)
pdf(file="~/Dropbox/conStruct/writeup/figs/bears/bear_std_xval.pdf",width=10,height=5,pointsize=14)
#	quartz(width=10,height=5)
	par(mfrow=c(1,2),mar=c(4,5,4,2))
	plot.xval.CIs(xval.CIs,K,jitter=0.1,xlim=c(0.75,7.25))
		legend(x="bottomright",pch=19,col=c("blue","green"),legend=c("spatial","nonspatial"))
	mtext("Predictive accuracy",side=2,padj=-5,font=2)
	plot.xval.CIs(xval.CIs,K,ylim=c(-200,0),jitter=0.1,xlim=c(3.75,7.25))
		legend(x="bottomright",pch=c(19,NA),lty=c(NA,1),lwd=c(NA,2),col=c(1,adjustcolor(1,0.8)),legend=c("mean","95% CI"))
	mtext("number of clusters",side=1,adj=-1.1,padj=4,font=2)
	mtext("Cross-validation results (Bears)",side=3,adj=5,padj=-2.5,font=2,cex=1.2)
dev.off()


library(conStruct)
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
load("~/Dropbox/gid_runs/mc_runs/bears/rep1_training.runs.sp.Robj")
for(k in 2:7){
	data.block$K <- k
	csr <- training.runs.sp[[k]][[1]]
	clst.match <- NULL
	if(k < 4){
		csr1.order <- NULL
	}
	clst.match$cols <- c(4,2)
	if(k > 2){
		tmp.csr <- training.runs.sp[[k-1]][[1]]
		clst.match <- match.clusters.x.runs(tmp.csr,csr,csr1.order)
		csr1.order <- clst.match$clst.order
	}
	pdf(file=paste0("~/Dropbox/conStruct/writeup/figs/bears/bear_sp",k,".pdf"),width=7,height=7)
		make.bear.redux.result.plot(csr = training.runs.sp[[k]][[1]],
									coords = bear.dataset$sample.coords,
									lump.dist = 200,
									cluster.colors = clst.match$cols,
									cluster.order=clst.match$clst.order,
									layout = matrix(c(rep(1,10),rep(2,15)),nrow=5,ncol=5,byrow=TRUE),
									box=TRUE)
	dev.off()
}

load("~/Dropbox/gid_runs/mc_runs/bears/rep1_training.runs.nsp.Robj")
for(k in 2:7){
	data.block$K <- k
	csr <- training.runs.nsp[[k]][[1]]
	clst.match <- NULL
	if(k < 4){
		csr1.order <- NULL
	}
	clst.match$cols <- c(4,2)
	if(k > 2){
		tmp.csr <- training.runs.nsp[[k-1]][[1]]
		clst.match <- match.clusters.x.runs(tmp.csr,csr,csr1.order)
		csr1.order <- clst.match$clst.order
	}
	pdf(file=paste0("~/Dropbox/conStruct/writeup/figs/bears/bear_nsp",k,".pdf"),width=7,height=7)
		make.bear.redux.result.plot(csr = training.runs.nsp[[k]][[1]],
									coords = bear.dataset$sample.coords,
									lump.dist = 200,
									cluster.colors = clst.match$cols,
									cluster.order=clst.match$clst.order,
									layout = matrix(c(rep(1,10),rep(2,15)),nrow=5,ncol=5,byrow=TRUE),
									box=TRUE)
	dev.off()
}

pdf(file=paste0("~/Dropbox/conStruct/writeup/figs/bears/bear_sp_clst_covs.pdf"),width=12,height=8,pointsize=14)
	#quartz(width=12,height=8)
	layout(matrix(c(1:6),nrow=2,ncol=3,byrow=TRUE))
	par(mar=c(4,5,3,2),oma=c(3,3,3,1))	
	plot.K.cluster.curves(K=7,data.block=data.block,training.runs=training.runs.sp,col.mat1=NULL,col.mat2=NULL)
	mtext(text="geographic distance",side=1,font=2,cex.axis=2,padj=4.5,adj=-3.55)
	mtext(text="allele frequency covariance",side=2,font=2,cex.axis=2,padj=-59.5,adj=20)
dev.off()


pdf(file=paste0("~/Dropbox/conStruct/writeup/figs/bears/bear_nsp_clst_covs.pdf"),width=12,height=8,pointsize=14)
	#quartz(width=12,height=8)
	layout(matrix(c(1:6),nrow=2,ncol=3,byrow=TRUE))
	par(mar=c(4,5,3,2),oma=c(3,3,3,1))	
	plot.K.cluster.curves(K=7,data.block=data.block,training.runs=training.runs.nsp,col.mat1=NULL,col.mat2=NULL)
	mtext(text="geographic distance",side=1,font=2,cex.axis=2,padj=4.5,adj=-3.55)
	mtext(text="allele frequency covariance",side=2,font=2,cex.axis=2,padj=-59.5,adj=20)
dev.off()


