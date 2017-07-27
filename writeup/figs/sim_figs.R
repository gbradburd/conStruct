################################################################
################################################################
#	Sim Figures for conStruct paper
################################################################
################################################################

source("~/Dropbox/conStruct/sims/cross_validation/summarize.xvals.R")

#K1
setwd("~/Dropbox/conStruct/sims/cross_validation/K_1/x_validation")
n.reps <- 10
K <- 7
for(n in 1:n.reps){
	load(sprintf("simK1_rep%s_test.lnl.Robj",n))
	assign(paste0("tl",n),test.lnl)
}
x.vals <- lapply(1:n.reps,function(n){get(sprintf("tl%s",n))})
x.vals.std <- lapply(x.vals,standardize.xvals)
xval.CIs <- get.xval.CIs(x.vals.std,K)
pdf(file="~/Dropbox/conStruct/writeup/figs/sims/simK1_std_xval.pdf",width=10,height=5,pointsize=14)
#	quartz(width=10,height=5)
	par(mfrow=c(1,2),mar=c(4,5,4,2))
	plot.xval.CIs(xval.CIs,K,simK=" (K = 1)")
		legend(x="bottomright",pch=19,col=c("blue","green"),legend=c("spatial","nonspatial"))
	mtext("Predictive accuracy",side=2,padj=-5)
	plot.xval.CIs(xval.CIs,K,simK=" (K = 1)",ylim=c(-10,0))
		legend(x="bottomleft",pch=c(19,NA),lty=c(NA,1),lwd=c(NA,2),col=c(1,adjustcolor(1,0.8)),legend=c("mean","95% CI"))
	mtext("number of clusters",side=1,adj=-0.95,padj=4)
	mtext("Cross-validation results (K=1)",side=3,adj=70,padj=-2.5,font=2,cex=1.2)
dev.off()


library(conStruct)
for(k in 2:7){
	load(sprintf("simK1__nsp_rep1K%s_data.block.Robj",k))
	load(sprintf("simK1__nsp_rep1K%s_conStruct.results.Robj",k))
	csr <- conStruct.results[[1]]
	clst.match <- NULL
	if(k < 4){
		csr1.order <- NULL
	}
	clst.match$cols <- c(4,2)
	if(k > 2){
		load(sprintf("simK1__nsp_rep1K%s_conStruct.results.Robj",k-1))
		tmp.csr <- conStruct.results[[1]]
		clst.match <- match.clusters.x.runs(tmp.csr,csr,csr1.order)
		csr1.order <- clst.match$clst.order
	}
	# pdf(file=paste0("~/Dropbox/conStruct/writeup/figs/sims/simK1_nsp_str_K",k,".pdf"),width=10,height=5)
		# make.structure.plot(data.block,csr,cluster.order=clst.match$clst.order,cluster.colors=clst.match$cols)
	# dev.off()
	pdf(file=paste0("~/Dropbox/conStruct/writeup/figs/sims/simK1_nsp_pies_K",k,".pdf"),width=5,height=5)
		make.admix.pie.plot(data.block,csr,cluster.colors=clst.match$cols,stat="MAP",title="",radii=3.5,x.lim=c(2.5,8.5),y.lim=c(2.5,8.5))
	dev.off()
}

library(conStruct)
for(k in 2:7){
	load(sprintf("simK1__sp_rep1K%s_data.block.Robj",k))
	load(sprintf("simK1__sp_rep1K%s_conStruct.results.Robj",k))
	csr <- conStruct.results[[1]]
	clst.match <- NULL
	if(k < 4){
		csr1.order <- NULL
	}
	clst.match$cols <- c(4,2)
	if(k > 2){
		load(sprintf("simK1__sp_rep1K%s_conStruct.results.Robj",k-1))
		tmp.csr <- conStruct.results[[1]]
		clst.match <- match.clusters.x.runs(tmp.csr,csr,csr1.order)
		csr1.order <- clst.match$clst.order
	}
	# pdf(file=paste0("~/Dropbox/conStruct/writeup/figs/sims/simK1_nsp_str_K",k,".pdf"),width=10,height=5)
		# make.structure.plot(data.block,csr,cluster.order=clst.match$clst.order,cluster.colors=clst.match$cols)
	# dev.off()
	pdf(file=paste0("~/Dropbox/conStruct/writeup/figs/sims/simK1_sp_pies_K",k,".pdf"),width=5,height=5)
		make.admix.pie.plot(data.block,csr,cluster.colors=clst.match$cols,stat="MAP",title="",radii=3.5,x.lim=c(2.5,8.5),y.lim=c(2.5,8.5))
	dev.off()
}
	pdf(file=paste0("~/Dropbox/conStruct/writeup/figs/sims/simK1_sp_pies_K",7,".pdf"),width=5,height=5)
		make.admix.pie.plot(data.block,csr,cluster.colors=c("blue","red","green","yellow","purple","orange","lightblue","darkgreen","lightblue","gray"),stat="MAP",title="",radii=3.5,x.lim=c(2.5,8.5),y.lim=c(2.5,8.5))
	dev.off()

#K2
setwd("~/Dropbox/conStruct/sims/cross_validation/K_2/xvals")
n.reps <- 10
K <- 7
for(n in 1:n.reps){
	load(sprintf("simK2_rep%s_test.lnl.Robj",n))
	assign(paste0("tl",n),test.lnl)
}
x.vals <- lapply(1:n.reps,function(n){get(sprintf("tl%s",n))})
x.vals.std <- lapply(x.vals,standardize.xvals)
xval.CIs <- get.xval.CIs(x.vals.std,K)
pdf(file="~/Dropbox/conStruct/writeup/figs/sims/simK2_std_xval.pdf",width=10,height=5,pointsize=14)
#	quartz(width=10,height=5)
	par(mfrow=c(1,2),mar=c(4,5,4,2))
	plot.xval.CIs(xval.CIs,K)
		mtext("Predictive accuracy",side=2,padj=-5)
	rect(1.7,-700,7.2,400,lty=2)
	arrows(1.7,-700,2.5,-3e3,lty=1,length=0.1)
	arrows(7.2,-700,6.6,-3e3,lty=1,length=0.1)
	TeachingDemos::subplot(fun = {
						plot.xval.CIs(xval.CIs,K,k.range=c(2:7),axes=FALSE,cex=1)
							axis(1,at=2:7,labels=c(2,"","","","",7),cex.axis=0.8,lty=2)
							axis(2,at=seq(-275,0,length.out=6),labels=c(-275,"","","","",0),cex.axis=0.8,lty=2)
							box(lwd=1.2,lty=2)
							legend(x="bottomright",pch=19,col=c("blue","green"),legend=c("spatial","nonspatial"),cex=0.8)
						},
					x=c(2.5,6.6),y=c(-12.5e3,-3e3))
		plot.xval.CIs(xval.CIs,K,simK=" (K = 2)",k.range=c(2:7),ylim=c(-9,0))
		legend(x="bottomleft",pch=c(19,NA),lty=c(NA,1),lwd=c(NA,2),col=c(1,adjustcolor(1,0.8)),legend=c("mean","95% CI"))
	mtext("number of clusters",side=1,adj=-0.95,padj=4)
	mtext("Cross-validation results (K=2)",side=3,adj=70,padj=-2.5,font=2,cex=1.2)
dev.off()

library(conStruct)
for(k in 2:7){
	load(sprintf("simK2__nsp_rep2K%s_data.block.Robj",k))
	load(sprintf("simK2__nsp_rep2K%s_conStruct.results.Robj",k))
	csr <- conStruct.results[[1]]
	clst.match <- NULL
	if(k < 4){
		csr1.order <- NULL
	}
	clst.match$cols <- c(4,2)
	if(k > 2){
		load(sprintf("simK2__nsp_rep2K%s_conStruct.results.Robj",k-1))
		tmp.csr <- conStruct.results[[1]]
		clst.match <- match.clusters.x.runs(tmp.csr,csr,csr1.order)
		csr1.order <- clst.match$clst.order
	}
	# pdf(file=paste0("~/Dropbox/conStruct/writeup/figs/sims/simK1_nsp_str_K",k,".pdf"),width=10,height=5)
		# make.structure.plot(data.block,csr,cluster.order=clst.match$clst.order,cluster.colors=clst.match$cols)
	# dev.off()
	pdf(file=paste0("~/Dropbox/conStruct/writeup/figs/sims/simK2_nsp_pies_K",k,".pdf"),width=5,height=5)
		make.admix.pie.plot(data.block,csr,cluster.colors=clst.match$cols,stat="MAP",title="",radii=3.5,x.lim=c(2.5,8.5),y.lim=c(2.5,8.5))
	dev.off()
}

library(conStruct) 
for(k in 2:7){
	load(sprintf("simK2__sp_rep2K%s_data.block.Robj",k))
	load(sprintf("simK2__sp_rep2K%s_conStruct.results.Robj",k))
	csr <- conStruct.results[[1]]
	clst.match <- NULL
	if(k < 4){
		csr1.order <- NULL
	}
	clst.match$cols <- c(4,2)
	if(k > 2){
		load(sprintf("simK2__sp_rep2K%s_conStruct.results.Robj",k-1))
		tmp.csr <- conStruct.results[[1]]
		clst.match <- match.clusters.x.runs(tmp.csr,csr,csr1.order)
		csr1.order <- clst.match$clst.order
	}
	# pdf(file=paste0("~/Dropbox/conStruct/writeup/figs/sims/simK1_nsp_str_K",k,".pdf"),width=10,height=5)
		# make.structure.plot(data.block,csr,cluster.order=clst.match$clst.order,cluster.colors=clst.match$cols)
	# dev.off()
	pdf(file=paste0("~/Dropbox/conStruct/writeup/figs/sims/simK2_sp_pies_K",k,".pdf"),width=5,height=5)
		make.admix.pie.plot(data.block,csr,cluster.colors=clst.match$cols,stat="MAP",title="",radii=3.5,x.lim=c(2.5,8.5),y.lim=c(2.5,8.5))
	dev.off()
}

setwd("~/Dropbox/conStruct/sims/cross_validation/K_3/xvals")
n.reps <- 10
K <- 7
for(n in 1:n.reps){
	load(sprintf("simK3_rep%s_test.lnl.Robj",n))
	assign(paste0("tl",n),test.lnl)
}
x.vals <- lapply(1:n.reps,function(n){get(sprintf("tl%s",n))})
x.vals.std <- lapply(x.vals,standardize.xvals)
xval.CIs <- get.xval.CIs(x.vals.std,K)

pdf(file="~/Dropbox/conStruct/writeup/figs/sims/simK3_std_xval.pdf",width=10,height=5,pointsize=14)
#	quartz(width=10,height=5)
	par(mfrow=c(1,2),mar=c(4,5,4,2))
	plot.xval.CIs(xval.CIs,K)
		mtext("Predictive accuracy",side=2,padj=-5)
	rect(2.8,-680,7.2,470,lty=2)
	arrows(2.8,-680,3.2,-4e3,lty=1,length=0.1)
	arrows(7.2,-700,6.6,-4e3,lty=1,length=0.1)
	TeachingDemos::subplot(fun = {
						plot.xval.CIs(xval.CIs,K,k.range=c(3:7),axes=FALSE,cex=1)
							axis(1,at=3:7,labels=c(3,"","","",7),cex.axis=0.8,lty=2)
							axis(2,at=seq(-140,0,length.out=6),labels=c(-140,"","","","",0),cex.axis=0.8,lty=2)
							box(lwd=1.2,lty=2)
							legend(x="bottomright",pch=19,col=c("blue","green"),legend=c("spatial","nonspatial"),cex=0.6)
						},
					x=c(3.2,6.6),y=c(-1.5e4,-4e3))
		plot.xval.CIs(xval.CIs,K,simK=" (K = 3)",k.range=c(3:7),ylim=c(-9,0))
		legend(x="bottomleft",pch=c(19,NA),lty=c(NA,1),lwd=c(NA,2),col=c(1,adjustcolor(1,0.8)),legend=c("mean","95% CI"))
	mtext("number of clusters",side=1,adj=-0.95,padj=4)
	mtext("Cross-validation results (K=3)",side=3,adj=70,padj=-2.5,font=2,cex=1.2)
dev.off()

library(conStruct)
load("~/Dropbox/conStruct/sims/cross_validation/K_3/sim.dataset.Robj")
for(k in 2:7){
	load(sprintf("simK3__nsp_rep2K%s_conStruct.results.Robj",k))
	csr <- conStruct.results[[1]]
	clst.match <- NULL
	if(k < 4){
		csr1.order <- NULL
	}
	clst.match$cols <- c(4,2)
	if(k > 2){
		load(sprintf("simK3__nsp_rep2K%s_conStruct.results.Robj",k-1))
		tmp.csr <- conStruct.results[[1]]
		clst.match <- match.clusters.x.runs(tmp.csr,csr,csr1.order)
		csr1.order <- clst.match$clst.order
	}
	# pdf(file=paste0("~/Dropbox/conStruct/writeup/figs/sims/simK1_nsp_str_K",k,".pdf"),width=10,height=5)
		# make.structure.plot(data.block,csr,cluster.order=clst.match$clst.order,cluster.colors=clst.match$cols)
	# dev.off()
	pdf(file=paste0("~/Dropbox/conStruct/writeup/figs/sims/simK3_nsp_pies_K",k,".pdf"),width=5,height=5)
		make.admix.pie.plot.tmp(sim.dataset$coords,k,sim.dataset$N,csr,cluster.colors=clst.match$cols,stat="MAP",title="",radii=3.5,x.lim=c(2.5,8.5),y.lim=c(2.5,8.5))
	dev.off()
}

library(conStruct)
load("~/Dropbox/conStruct/sims/cross_validation/K_3/sim.dataset.Robj")
for(k in 2:7){
	load(sprintf("simK3__sp_rep2K%s_conStruct.results.Robj",k))
	csr <- conStruct.results[[1]]
	clst.match <- NULL
	if(k < 4){
		csr1.order <- NULL
	}
	clst.match$cols <- c(4,2)
	if(k > 2){
		load(sprintf("simK3__sp_rep2K%s_conStruct.results.Robj",k-1))
		tmp.csr <- conStruct.results[[1]]
		clst.match <- match.clusters.x.runs(tmp.csr,csr,csr1.order)
		csr1.order <- clst.match$clst.order
	}
	# pdf(file=paste0("~/Dropbox/conStruct/writeup/figs/sims/simK1_nsp_str_K",k,".pdf"),width=10,height=5)
		# make.structure.plot(data.block,csr,cluster.order=clst.match$clst.order,cluster.colors=clst.match$cols)
	# dev.off()
	pdf(file=paste0("~/Dropbox/conStruct/writeup/figs/sims/simK3_sp_pies_K",k,".pdf"),width=5,height=5)
		make.admix.pie.plot.tmp(sim.dataset$coords,k,sim.dataset$N,csr,cluster.colors=clst.match$cols,stat="MAP",title="",radii=3.5,x.lim=c(2.5,8.5),y.lim=c(2.5,8.5))
	dev.off()
}





#All together now

#K1
setwd("~/Dropbox/conStruct/sims/cross_validation/K_1/x_validation")
n.reps <- 10
K <- 7
for(n in 1:n.reps){
	load(sprintf("simK1_rep%s_test.lnl.Robj",n))
	assign(paste0("tl",n),test.lnl)
}
x.vals <- lapply(1:n.reps,function(n){get(sprintf("tl%s",n))})
x.vals.std <- lapply(x.vals,standardize.xvals)
xval.CIs <- get.xval.CIs(x.vals.std,K)
pdf(file="~/Dropbox/conStruct/writeup/figs/sims/xvals.pdf",width=10,height=15,pointsize=14)
#	quartz(width=10,height=15)
	layout(matrix(c(1:6),nrow=3,ncol=2,byrow=TRUE))
	par(mar=c(4,5,3,2),oma=c(3,3,3,1))
	plot.xval.CIs(xval.CIs,K,simK=" (K = 1)")
		legend(x="bottomright",pch=19,col=c("blue","green"),legend=c("spatial","nonspatial"))
	plot.xval.CIs(xval.CIs,K,simK=" (K = 1)",ylim=c(-10,0))
		legend(x="bottomleft",pch=c(19,NA),lty=c(NA,1),lwd=c(NA,2),col=c(1,adjustcolor(1,0.8)),legend=c("mean","95% CI"))
	mtext("K = 1",side=3,adj=-0.3,padj=-2.2,font=2,cex=1.2)


#K2
setwd("~/Dropbox/conStruct/sims/cross_validation/K_2/xvals")
n.reps <- 10
K <- 7
for(n in 1:n.reps){
	load(sprintf("simK2_rep%s_test.lnl.Robj",n))
	assign(paste0("tl",n),test.lnl)
}
x.vals <- lapply(1:n.reps,function(n){get(sprintf("tl%s",n))})
x.vals.std <- lapply(x.vals,standardize.xvals)
xval.CIs <- get.xval.CIs(x.vals.std,K)
#	quartz(width=10,height=5)
	plot.xval.CIs(xval.CIs,K)
		mtext("Predictive accuracy",side=2,padj=-5,font=2)
	rect(1.7,-700,7.2,400,lty=2)
	arrows(1.7,-700,2.5,-3e3,lty=1,length=0.1)
	arrows(7.2,-700,6.6,-3e3,lty=1,length=0.1)
	TeachingDemos::subplot(fun = {
						plot.xval.CIs(xval.CIs,K,k.range=c(2:7),axes=FALSE,cex=1)
							axis(1,at=2:7,labels=c(2,"","","","",7),cex.axis=0.8,lty=2)
							axis(2,at=seq(-275,0,length.out=6),labels=c(-275,"","","","",0),cex.axis=0.8,lty=2)
							box(lwd=1.2,lty=2)
						},
					x=c(2.5,6.6),y=c(-12.5e3,-3e3))
		plot.xval.CIs(xval.CIs,K,simK=" (K = 2)",k.range=c(2:7),ylim=c(-9,0))
	mtext("K = 2",side=3,adj=-0.3,padj=-2.2,font=2,cex=1.2)
	
setwd("~/Dropbox/conStruct/sims/cross_validation/K_3/xvals")
n.reps <- 10
K <- 7
for(n in 1:n.reps){
	load(sprintf("simK3_rep%s_test.lnl.Robj",n))
	assign(paste0("tl",n),test.lnl)
}
x.vals <- lapply(1:n.reps,function(n){get(sprintf("tl%s",n))})
x.vals.std <- lapply(x.vals,standardize.xvals)
xval.CIs <- get.xval.CIs(x.vals.std,K)

	plot.xval.CIs(xval.CIs,K)
	rect(2.8,-680,7.2,470,lty=2)
	arrows(2.8,-680,3.2,-4e3,lty=1,length=0.1)
	arrows(7.2,-700,6.6,-4e3,lty=1,length=0.1)
	TeachingDemos::subplot(fun = {
						plot.xval.CIs(xval.CIs,K,k.range=c(3:7),axes=FALSE,cex=1)
							axis(1,at=3:7,labels=c(3,"","","",7),cex.axis=0.8,lty=2)
							axis(2,at=seq(-140,0,length.out=6),labels=c(-140,"","","","",0),cex.axis=0.8,lty=2)
							box(lwd=1.2,lty=2)
						},
					x=c(3.2,6.6),y=c(-1.5e4,-4e3))
		plot.xval.CIs(xval.CIs,K,simK=" (K = 3)",k.range=c(3:7),ylim=c(-9,0))
	mtext("K = 3",side=3,adj=-0.3,padj=-2.2,font=2,cex=1.2)
	mtext("number of clusters",side=1,adj=-0.7,padj=4,font=2)
dev.off()

