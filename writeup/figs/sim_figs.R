################################################################
################################################################
#	Sim Figures for conStruct paper
################################################################
################################################################
source("~/Dropbox/conStruct/writeup/figs/fig_funcs.R")

#K1
setwd("~/Dropbox/conStruct/sims/cross_validation/K_1/x_validation")
pdf(file="~/Dropbox/conStruct/writeup/figs/sims/simK1_std_xval.pdf",width=10,height=5,pointsize=14)
#	quartz(width=10,height=5)
	par(mfrow=c(1,2),mar=c(4,5,4,2))
	plot.sim.xvals(dir="~/Dropbox/conStruct/sims/cross_validation/K_1/x_validation",n.reps=10,K=7,simK=1,y.lim=c(-10,0))
dev.off()

setwd("~/Dropbox/conStruct/sims/cross_validation/K_1")
library(conStruct)
load("sim.dataset.Robj")
freq.data <- conStruct:::process.freq.data(sim.dataset$freq.data$freqs)
data.block <- conStruct:::make.data.block(K = 1,
										  freq.data = freq.data,
										  coords = sim.dataset$coords,
										  spatial = TRUE,
										  geoDist = fields::rdist(sim.dataset$coords))

output.list.sp <- vector("list",7)
for(k in 1:7){
	load(sprintf("~/Dropbox/conStruct/sims/cross_validation/K_1/runs/simK1_K%s_sp_conStruct.results.Robj",k))
	output.list.sp[[k]] <- conStruct.results
}

output.list.nsp <- vector("list",7)
for(k in 1:7){
	load(sprintf("~/Dropbox/conStruct/sims/cross_validation/K_1/runs/simK1_K%s_nsp_conStruct.results.Robj",k))
	output.list.nsp[[k]] <- conStruct.results
}


plot.sim.pies(data.block = data.block,
			  K = 7,
			  output.list = output.list.nsp,
			  file.name = "~/Dropbox/conStruct/writeup/figs/sims/simK1_nsp_pies_K")

plot.sim.pies(data.block = data.block,
			  K = 7,
			  output.list = output.list.sp,
			  file.name = "~/Dropbox/conStruct/writeup/figs/sims/simK1_sp_pies_K")

K <- 7
laycon.sp <- matrix(0,nrow=7,ncol=7)
colnames(laycon.sp) <- paste(1:7)
for(k in 1:K){
	csr <- output.list.sp[[k]][[1]]
	clst.match <- NULL
	if(k < 4){
		csr1.order <- NULL
	}
	if(k > 2){
		tmp.csr <- output.list.sp[[k-1]][[1]]
		clst.match <- match.clusters.x.runs(tmp.csr,csr,csr1.order)
		csr1.order <- clst.match$clst.order
	}
	if(is.null(csr1.order)){
		csr1.order <- 1:k
	}
	laycon.sp[1:k,k] <- calculate.layer.contributions(csr$MAP,
														data.block)[csr1.order]
}

laycon.nsp <- matrix(0,nrow=7,ncol=7)
colnames(laycon.nsp) <- paste(1:7)
for(k in 1:K){
	csr <- output.list.nsp[[k]][[1]]
	clst.match <- NULL
	if(k <= 2){
		csr1.order <- NULL
	}
	if(k==2){
		tmp.csr <- output.list.sp[[k]][[1]]
		clst.match <- match.clusters.x.runs(tmp.csr,csr,csr1.order)
		csr1.order <- clst.match$clst.order
	}
	if(k > 2){
		tmp.csr <- output.list.nsp[[k-1]][[1]]
		clst.match <- match.clusters.x.runs(tmp.csr,csr,csr1.order)
		csr1.order <- clst.match$clst.order
	}
	if(is.null(csr1.order)){
		csr1.order <- 1:k
	}
	laycon.nsp[1:k,k] <- calculate.layer.contributions(csr$MAP,
														data.block)[csr1.order]
}

pdf(file="~/Dropbox/conStruct/writeup/figs/sims/simK1_laycon_barplots.pdf",width=8,height=4,pointsize=14)
	par(mfrow=c(1,2),mar=c(4,4,3,0.5))
	barplot(laycon.sp,	
			col=c("blue","red","green","yellow","purple","orange","lightblue"),
			xlab="",ylab="layer importance")
	barplot(laycon.nsp,
			col=c("blue","red","green","yellow","purple","orange","lightblue"),
			xlab="",ylab="")
			mtext(side=1,text="number of layers",padj=4.25,adj=-1)
			mtext(side=3,text="Layer contributions (K=1)",padj=-2.25,adj=-21.5,font=2,cex=1.2)
dev.off()


#K2
setwd("~/Dropbox/conStruct/sims/cross_validation/K_2/x_validation")
pdf(file="~/Dropbox/conStruct/writeup/figs/sims/simK2_std_xval.pdf",width=10,height=5,pointsize=14)
#	quartz(width=10,height=5)
	par(mfrow=c(1,2),mar=c(4,5,4,2))
	plot.sim.xvals(dir="~/Dropbox/conStruct/sims/cross_validation/K_2/x_validation",n.reps=10,K=7,simK=2,y.lim=c(-10,0))
dev.off()

setwd("~/Dropbox/conStruct/sims/cross_validation/K_2")
library(conStruct)
load("sim.dataset.Robj")
freq.data <- conStruct:::process.freq.data(sim.dataset$freq.data$freqs)
data.block <- conStruct:::make.data.block(K = 1,
										  freq.data = freq.data,
										  coords = sim.dataset$coords,
										  spatial = TRUE,
										  geoDist = fields::rdist(sim.dataset$coords))

output.list.sp <- vector("list",7)
for(k in 1:7){
	load(sprintf("~/Dropbox/conStruct/sims/cross_validation/K_2/runs/simK2_K%s_sp_conStruct.results.Robj",k))
	output.list.sp[[k]] <- conStruct.results
}

output.list.nsp <- vector("list",7)
for(k in 1:7){
	load(sprintf("~/Dropbox/conStruct/sims/cross_validation/K_2/runs/simK2_K%s_nsp_conStruct.results.Robj",k))
	output.list.nsp[[k]] <- conStruct.results
}


plot.sim.pies(data.block = data.block,
			  K = 7,
			  output.list = output.list.nsp,
			  file.name = "~/Dropbox/conStruct/writeup/figs/sims/simK2_nsp_pies_K")

plot.sim.pies(data.block = data.block,
			  K = 7,
			  output.list = output.list.sp,
			  file.name = "~/Dropbox/conStruct/writeup/figs/sims/simK2_sp_pies_K")			

K <- 7
laycon.sp <- matrix(0,nrow=7,ncol=7)
colnames(laycon.sp) <- paste(1:7)
for(k in 1:K){
	csr <- output.list.sp[[k]][[1]]
	clst.match <- NULL
	if(k < 4){
		csr1.order <- NULL
	}
	if(k > 2){
		tmp.csr <- output.list.sp[[k-1]][[1]]
		clst.match <- match.clusters.x.runs(tmp.csr,csr,csr1.order)
		csr1.order <- clst.match$clst.order
	}
	if(is.null(csr1.order)){
		csr1.order <- 1:k
	}
	laycon.sp[1:k,k] <- calculate.layer.contributions(csr$MAP,
														data.block)[csr1.order]
}

laycon.nsp <- matrix(0,nrow=7,ncol=7)
colnames(laycon.nsp) <- paste(1:7)
for(k in 1:K){
	csr <- output.list.nsp[[k]][[1]]
	clst.match <- NULL
	if(k <= 2){
		csr1.order <- NULL
	}
	if(k==2){
		tmp.csr <- output.list.sp[[k]][[1]]
		clst.match <- match.clusters.x.runs(tmp.csr,csr,csr1.order)
		csr1.order <- clst.match$clst.order
	}
	if(k > 2){
		tmp.csr <- output.list.nsp[[k-1]][[1]]
		clst.match <- match.clusters.x.runs(tmp.csr,csr,csr1.order)
		csr1.order <- clst.match$clst.order
	}
	if(is.null(csr1.order)){
		csr1.order <- 1:k
	}
	laycon.nsp[1:k,k] <- calculate.layer.contributions(csr$MAP,
														data.block)[csr1.order]
}

pdf(file="~/Dropbox/conStruct/writeup/figs/sims/simK2_laycon_barplots.pdf",width=8,height=4,pointsize=14)
	par(mfrow=c(1,2),mar=c(4,4,3,0.5))
	barplot(laycon.sp,	
			col=c("blue","red","green","yellow","purple","orange","lightblue"),
			xlab="",ylab="layer importance")
	barplot(laycon.nsp,
			col=c("blue","red","green","yellow","purple","orange","lightblue"),
			xlab="",ylab="")
			mtext(side=1,text="number of layers",padj=4.25,adj=-1)
			mtext(side=3,text="Layer contributions (K=2)",padj=-2.25,adj=-21.5,font=2,cex=1.2)
dev.off()

pdf(file="~/Dropbox/conStruct/writeup/figs/sims/simK2_adprop_fit.pdf",width=6,height=6,pointsize=14)
	viz.admix.results(sim.admix.props = sim.dataset$admix.list$w,
				 	 conStruct.results = output.list.sp[[2]][[1]],
				 	 clst.order=c(2,1))
dev.off()

#K3
setwd("~/Dropbox/conStruct/sims/cross_validation/K_3/x_validation")
pdf(file="~/Dropbox/conStruct/writeup/figs/sims/simK3_std_xval.pdf",width=10,height=5,pointsize=14)
#	quartz(width=10,height=5)
	par(mfrow=c(1,2),mar=c(4,5,4,2))
	plot.sim.xvals(dir="~/Dropbox/conStruct/sims/cross_validation/K_3/x_validation",n.reps=10,K=7,simK=3,y.lim=c(-10,0))
dev.off()

setwd("~/Dropbox/conStruct/sims/cross_validation/K_3")
library(conStruct)
load("sim.dataset.Robj")
freq.data <- conStruct:::process.freq.data(sim.dataset$freq.data$freqs)
data.block <- conStruct:::make.data.block(K = 1,
										  freq.data = freq.data,
										  coords = sim.dataset$coords,
										  spatial = TRUE,
										  geoDist = fields::rdist(sim.dataset$coords))

output.list.sp <- vector("list",7)
for(k in 1:7){
	load(sprintf("~/Dropbox/conStruct/sims/cross_validation/K_3/runs/simK3_K%s_sp_conStruct.results.Robj",k))
	output.list.sp[[k]] <- conStruct.results
}

output.list.nsp <- vector("list",7)
for(k in 1:7){
	load(sprintf("~/Dropbox/conStruct/sims/cross_validation/K_3/runs/simK3_K%s_nsp_conStruct.results.Robj",k))
	output.list.nsp[[k]] <- conStruct.results
}


plot.sim.pies(data.block = data.block,
			  K = 7,
			  output.list = output.list.nsp,
			  file.name = "~/Dropbox/conStruct/writeup/figs/sims/simK3_nsp_pies_K")

plot.sim.pies(data.block = data.block,
			  K = 7,
			  output.list = output.list.sp,
			  file.name = "~/Dropbox/conStruct/writeup/figs/sims/simK3_sp_pies_K")		

K <- 7
laycon.sp <- matrix(0,nrow=7,ncol=7)
colnames(laycon.sp) <- paste(1:7)
for(k in 1:K){
	csr <- output.list.sp[[k]][[1]]
	clst.match <- NULL
	if(k < 4){
		csr1.order <- NULL
	}
	if(k > 2){
		tmp.csr <- output.list.sp[[k-1]][[1]]
		clst.match <- match.clusters.x.runs(tmp.csr,csr,csr1.order)
		csr1.order <- clst.match$clst.order
	}
	if(is.null(csr1.order)){
		csr1.order <- 1:k
	}
	laycon.sp[1:k,k] <- calculate.layer.contributions(csr$MAP,
														data.block)[csr1.order]
}

laycon.nsp <- matrix(0,nrow=7,ncol=7)
colnames(laycon.nsp) <- paste(1:7)
for(k in 1:K){
	csr <- output.list.nsp[[k]][[1]]
	clst.match <- NULL
	if(k <= 2){
		csr1.order <- NULL
	}
	if(k==2){
		tmp.csr <- output.list.sp[[k]][[1]]
		clst.match <- match.clusters.x.runs(tmp.csr,csr,csr1.order)
		csr1.order <- clst.match$clst.order
	}
	if(k > 2){
		tmp.csr <- output.list.nsp[[k-1]][[1]]
		clst.match <- match.clusters.x.runs(tmp.csr,csr,csr1.order)
		csr1.order <- clst.match$clst.order
	}
	if(is.null(csr1.order)){
		csr1.order <- 1:k
	}
	laycon.nsp[1:k,k] <- calculate.layer.contributions(csr$MAP,
														data.block)[csr1.order]
}

pdf(file="~/Dropbox/conStruct/writeup/figs/sims/simK3_laycon_barplots.pdf",width=8,height=4,pointsize=14)
	par(mfrow=c(1,2),mar=c(4,4,3,0.5))
	barplot(laycon.sp,	
			col=c("blue","red","green","yellow","purple","orange","lightblue"),
			xlab="",ylab="layer importance")
	barplot(laycon.nsp,
			col=c("blue","red","green","yellow","purple","orange","lightblue"),
			xlab="",ylab="")
			mtext(side=1,text="number of layers",padj=4.25,adj=-1)
			mtext(side=3,text="Layer contributions (K=3)",padj=-2.25,adj=-21.5,font=2,cex=1.2)
dev.off()

pdf(file="~/Dropbox/conStruct/writeup/figs/sims/simK3_adprop_fit.pdf",width=6,height=6,pointsize=14)
viz.admix.results(sim.admix.props = sim.dataset$admix.list$w,
				  conStruct.results = output.list.sp[[3]][[1]],
				  clst.order=c(3,2,1))
dev.off()


pdf(file="~/Dropbox/conStruct/writeup/figs/sims/sim_xvals.pdf",width=16,height=5,pointsize=14)
#	quartz(width=15,height=5)
#K1
setwd("~/Dropbox/gid_runs/mc_runs/sims/simK1")
n.reps <- 10
K <- 7
for(n in 1:n.reps){
	load(sprintf("simK1_rep%s_test.lnl.Robj",n))
	assign(paste0("tl",n),test.lnl)
}
x.vals <- lapply(1:n.reps,function(n){get(sprintf("tl%s",n))})
x.vals.std <- lapply(x.vals,standardize.xvals)
xval.CIs <- get.xval.CIs(x.vals.std,K)
	par(mfrow=c(1,3),mar=c(4.5,5,4,2))
	plot.xval.CIs(xval.CIs,K,xlim=c(0.8,7.2),jitter=0.15)
		mtext("Predictive accuracy",side=2,padj=-4)
	# rect(0.7,-120,7.2,20,lty=2)
	# arrows(0.7,-120,3,-3e2,lty=1,length=0.2)
	# arrows(7.2,-120,7,-3e2,lty=1,length=0.2)
	# TeachingDemos::subplot(fun = {
						# plot.xval.CIs(xval.CIs,K,k.range=c(1:7),axes=FALSE,cex=1,ylim=c(-120,0))
							# axis(1,at=1:7,labels=c(1,"","","","","",7),cex.axis=0.8,lty=2)
							# axis(2,at=seq(-120,0,length.out=6),labels=c(-120,"","","","",0),cex.axis=0.8,lty=2)
							# box(lwd=1.2,lty=2)
							# legend(x="bottomright",pch=19,col=c("blue","green"),legend=c("spatial","nonspatial"),cex=0.8)
						# },
					# x=c(3,7),y=c(-620,-3e2))
	legend(x="bottomright",pch=19,col=c("blue","green"),legend=c("spatial","nonspatial"),cex=1.2)
	mtext("K = 1",side=3,adj=0.5,padj=-1.5,font=2,cex=1.2)
#K2
setwd("~/Dropbox/gid_runs/mc_runs/sims/simK2")
n.reps <- 10
K <- 7
for(n in 1:n.reps){
	load(sprintf("simK2_rep%s_test.lnl.Robj",n))
	assign(paste0("tl",n),test.lnl)
}
x.vals <- lapply(1:n.reps,function(n){get(sprintf("tl%s",n))})
x.vals.std <- lapply(x.vals,standardize.xvals)
xval.CIs <- get.xval.CIs(x.vals.std,K)
	plot.xval.CIs(xval.CIs,K,xlim=c(0.8,7.2),jitter=0.15)
	rect(0.7,-700,7.3,400,lty=2)
	arrows(0.7,-700,2.5,-3e3,lty=1,length=0.15)
	arrows(7.2,-700,6.6,-3e3,lty=1,length=0.15)
	TeachingDemos::subplot(fun = {
						plot.xval.CIs(xval.CIs,K,k.range=c(1:7),xlim=c(0.7,7.2),ylim=c(-275,0),axes=FALSE,cex=1,jitter=0.1)
							axis(1,at=1:7,labels=c(1,"","","","","",7),cex.axis=0.8,lty=2)
							axis(2,at=seq(-275,0,length.out=6),labels=c(-275,"","","","",0),cex.axis=0.8,lty=2)
							box(lwd=1.2,lty=2)
						},
					x=c(2.5,6.6),y=c(-12.5e3,-3e3))
	mtext("K = 2",side=3,adj=0.5,padj=-1.5,font=2,cex=1.2)
	mtext("number of clusters",side=1,padj=3.4)
#K3
setwd("~/Dropbox/gid_runs/mc_runs/sims/simK3")
n.reps <- 10
K <- 7
for(n in 1:n.reps){
	load(sprintf("simK3_rep%s_test.lnl.Robj",n))
	assign(paste0("tl",n),test.lnl)
}
x.vals <- lapply(1:n.reps,function(n){get(sprintf("tl%s",n))})
x.vals.std <- lapply(x.vals,standardize.xvals)
xval.CIs <- get.xval.CIs(x.vals.std,K)
	plot.xval.CIs(xval.CIs,K,xlim=c(0.8,7.2),jitter=0.15)
	rect(0.7,-680,7.3,470,lty=2)
	arrows(0.7,-680,3.2,-4e3,lty=1,length=0.1)
	arrows(7.2,-700,6.6,-4e3,lty=1,length=0.1)
	TeachingDemos::subplot(fun = {
						plot.xval.CIs(xval.CIs,K,k.range=c(1:7),ylim=c(-140,0),xlim=c(0.7,7.3),axes=FALSE,cex=1,jitter=0.1)
							axis(1,at=1:7,labels=c(1,"","","","","",7),cex.axis=0.8,lty=2)
							axis(2,at=seq(-140,0,length.out=6),labels=c(-140,"","","","",0),cex.axis=0.8,lty=2)
							box(lwd=1.2,lty=2)
						},
					x=c(3.2,6.6),y=c(-1.5e4,-4e3))
	mtext("K = 3",side=3,adj=0.5,padj=-1.5,font=2,cex=1.2)
dev.off()

# #All together now

# #K1
# setwd("~/Dropbox/conStruct/sims/cross_validation/K_1/x_validation")
# n.reps <- 10
# K <- 7
# for(n in 1:n.reps){
	# load(sprintf("simK1_rep%s_test.lnl.Robj",n))
	# assign(paste0("tl",n),test.lnl)
# }
# x.vals <- lapply(1:n.reps,function(n){get(sprintf("tl%s",n))})
# x.vals.std <- lapply(x.vals,standardize.xvals)
# xval.CIs <- get.xval.CIs(x.vals.std,K)
# pdf(file="~/Dropbox/conStruct/writeup/figs/sims/xvals.pdf",width=10,height=15,pointsize=14)
# #	quartz(width=10,height=15)
	# layout(matrix(c(1:6),nrow=3,ncol=2,byrow=TRUE))
	# par(mar=c(4,5,3,2),oma=c(3,3,3,1))
	# plot.xval.CIs(xval.CIs,K,xlim=c(0.7,7.3))
		# legend(x="bottomright",pch=19,col=c("blue","green"),legend=c("spatial","nonspatial"))
	# rect(0.7,-110,7.3,10,lwd=1,lty=2)
	# plot.xval.CIs(xval.CIs,K,ylim=c(-100,0),xlim=c(0.7,7.3))
		# legend(x="bottomleft",pch=c(19,NA),lty=c(NA,1),lwd=c(NA,2),col=c(1,adjustcolor(1,0.8)),legend=c("mean","95% CI"))
	# mtext("K = 1",side=3,adj=-0.3,padj=-2.2,font=2,cex=1.2)

# #K2
# setwd("~/Dropbox/conStruct/sims/cross_validation/K_2/xvals")
# n.reps <- 10
# K <- 7
# for(n in 1:n.reps){
	# load(sprintf("simK2_rep%s_test.lnl.Robj",n))
	# assign(paste0("tl",n),test.lnl)
# }
# x.vals <- lapply(1:n.reps,function(n){get(sprintf("tl%s",n))})
# x.vals.std <- lapply(x.vals,standardize.xvals)
# xval.CIs <- get.xval.CIs(x.vals.std,K)
# #	quartz(width=10,height=5)
	# plot.xval.CIs(xval.CIs,K)
		# mtext("Predictive accuracy",side=2,padj=-5,font=2)
	# rect(1.7,-700,7.2,400,lty=2)
	# arrows(1.7,-700,2.5,-3e3,lty=1,length=0.1)
	# arrows(7.2,-700,6.6,-3e3,lty=1,length=0.1)
	# TeachingDemos::subplot(fun = {
						# plot.xval.CIs(xval.CIs,K,k.range=c(2:7),axes=FALSE,cex=1)
							# axis(1,at=2:7,labels=c(2,"","","","",7),cex.axis=0.8,lty=2)
							# axis(2,at=seq(-275,0,length.out=6),labels=c(-275,"","","","",0),cex.axis=0.8,lty=2)
							# box(lwd=1.2,lty=2)
						# },
					# x=c(2.5,6.6),y=c(-12.5e3,-3e3))
		# plot.xval.CIs(xval.CIs,K,simK=" (K = 2)",k.range=c(2:7),ylim=c(-9,0))
	# mtext("K = 2",side=3,adj=-0.3,padj=-2.2,font=2,cex=1.2)
	
# setwd("~/Dropbox/conStruct/sims/cross_validation/K_3/xvals")
# n.reps <- 10
# K <- 7
# for(n in 1:n.reps){
	# load(sprintf("simK3_rep%s_test.lnl.Robj",n))
	# assign(paste0("tl",n),test.lnl)
# }
# x.vals <- lapply(1:n.reps,function(n){get(sprintf("tl%s",n))})
# x.vals.std <- lapply(x.vals,standardize.xvals)
# xval.CIs <- get.xval.CIs(x.vals.std,K)

	# plot.xval.CIs(xval.CIs,K)
	# rect(2.8,-680,7.2,470,lty=2)
	# arrows(2.8,-680,3.2,-4e3,lty=1,length=0.1)
	# arrows(7.2,-700,6.6,-4e3,lty=1,length=0.1)
	# TeachingDemos::subplot(fun = {
						# plot.xval.CIs(xval.CIs,K,k.range=c(3:7),axes=FALSE,cex=1)
							# axis(1,at=3:7,labels=c(3,"","","",7),cex.axis=0.8,lty=2)
							# axis(2,at=seq(-140,0,length.out=6),labels=c(-140,"","","","",0),cex.axis=0.8,lty=2)
							# box(lwd=1.2,lty=2)
						# },
					# x=c(3.2,6.6),y=c(-1.5e4,-4e3))
		# plot.xval.CIs(xval.CIs,K,simK=" (K = 3)",k.range=c(3:7),ylim=c(-9,0))
	# mtext("K = 3",side=3,adj=-0.3,padj=-2.2,font=2,cex=1.2)
	# mtext("number of clusters",side=1,adj=-0.7,padj=4,font=2)
# dev.off()

