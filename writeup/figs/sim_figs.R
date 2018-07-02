################################################################
################################################################
#	Sim Figures for conStruct paper
################################################################
################################################################

source("fig_funcs.R")

################################
#K1
################################
pdf(file="sims/simK1_std_xval.pdf",width=10,height=5,pointsize=14)
#	quartz(width=10,height=5)
	par(mfrow=c(1,2),mar=c(4,5,4,2))
	plot.sim.xvals(dir="../../sims/analyses/conStruct/cross_validation/K_1",n.reps=10,K=7,simK=1,y.lim=c(-10,0))
dev.off()


load("../../sims/sim_data/K_1/sim.dataset.Robj")
freq.data <- conStruct:::process.freq.data(sim.dataset$freq.data$freqs)
data.block <- conStruct:::make.data.block(K = 1,
										  freq.data = freq.data,
										  coords = sim.dataset$coords,
										  spatial = TRUE,
										  geoDist = fields::rdist(sim.dataset$coords))

output.list.sp <- vector("list",7)
for(k in 1:7){
	load(sprintf("../../sims/analyses/conStruct/runs/K_1/simK1_K%s_sp_conStruct.results.Robj",k))
	conStruct.results <- cluster.2.layer(conStruct.results)
	for(j in 1:k){
		names(conStruct.results[[1]]$MAP$layer.params[[j]])[4] <- "phi"
	}
	output.list.sp[[k]] <- conStruct.results
}

	output.list.nsp <- vector("list",7)
	for(k in 1:7){
		load(sprintf("../../sims/analyses/conStruct/runs/K_1/simK1_K%s_nsp_conStruct.results.Robj",k))
		conStruct.results <- cluster.2.layer(conStruct.results)
		for(j in 1:k){
			names(conStruct.results[[1]]$MAP$layer.params[[j]])[4] <- "phi"
		}
		output.list.nsp[[k]] <- conStruct.results
	}

pdf(file="sims/simK1_nsp_pies.pdf",width=8,height=6.3,pointsize=14)
	par(mfrow=c(2,3),oma=c(1,0,3.5,0))
	plot.sim.pies.multipanel(data.block = data.block,
							 K = 7,
							 output.list = output.list.nsp,
							 radii = 1.7,
							 mar = c(4,2,1,2),trueK=1)
dev.off()

pdf(file="sims/simK1_sp_pies.pdf",width=8,height=6.3,pointsize=14)
	par(mfrow=c(2,3),oma=c(1,0,3.5,0))
	plot.sim.pies.multipanel(data.block = data.block,
							 K = 7,
							 output.list = output.list.sp,
							 radii = 1.7,
							 mar = c(4,2,1,2),trueK=1)
dev.off()

pdf(file="sims/Fig2_simK1_sp_vs_nsp.pdf",width=8,height=6.3,pointsize=14)
	#quartz(width=6,height=4,pointsize=14)
	radii <- 1.7
	mar <- c(4,2,1,2)
	par(mfrow=c(2,3),oma=c(1,0,3.5,0))
		make.admix.pie.plot(output.list.nsp[[2]][[1]]$MAP$admix.proportions,
								data.block$coords,radii=radii,x.lim=c(2.5,8.5),y.lim=c(2.5,8.5),mar= mar)
			mtext(side=1,text=expression(paste("(a) ",italic("K")," = 2")),padj=2.7,adj=0.4)
		make.admix.pie.plot(output.list.nsp[[3]][[1]]$MAP$admix.proportions[,c(2,3,1)],
								data.block$coords,radii=radii,x.lim=c(2.5,8.5),y.lim=c(2.5,8.5),mar= mar)
			mtext(side=1,text=expression(paste("(b) ",italic("K")," = 3")),padj=2.7,adj=0.4)
			mtext(side=3,text=bquote(paste("True ",italic("K")," = 1")),cex=1.5,padj=-0.7)
		make.admix.pie.plot(output.list.nsp[[4]][[1]]$MAP$admix.proportions[,c(4,2,3,1)],
								data.block$coords,radii=radii,x.lim=c(2.5,8.5),y.lim=c(2.5,8.5),mar= mar)
			mtext(side=1,text=expression(paste("(c) ",italic("K")," = 4")),padj=2.7,adj=0.4)
		make.admix.pie.plot(output.list.sp[[2]][[1]]$MAP$admix.proportions,
								data.block$coords,radii=radii,x.lim=c(2.5,8.5),y.lim=c(2.5,8.5),mar= mar)
			mtext(side=1,text=expression(paste("(d) ",italic("K")," = 2")),padj=2.7,adj=0.4)
		make.admix.pie.plot(output.list.sp[[3]][[1]]$MAP$admix.proportions[,c(3,1,2)],
								data.block$coords,radii=radii,x.lim=c(2.5,8.5),y.lim=c(2.5,8.5),mar= mar)
			mtext(side=1,text=expression(paste("(e) ",italic("K")," = 3")),padj=2.7,adj=0.4)
		make.admix.pie.plot(output.list.sp[[4]][[1]]$MAP$admix.proportions[,c(4,2,1,3)],
								data.block$coords,radii=radii,x.lim=c(2.5,8.5),y.lim=c(2.5,8.5),mar= mar)
			mtext(side=1,text=expression(paste("(f) ",italic("K")," = 4")),padj=2.7,adj=0.4)
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
		csr1.order <- match.layers.x.runs(tmp.csr$MAP$admix.proportions,csr$MAP$admix.proportions,csr1.order)
	}
	if(is.null(csr1.order)){
		csr1.order <- 1:k
	}
	laycon.sp[1:k,k] <- calculate.layer.contribution(csr,data.block,csr1.order)
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
		csr1.order <- match.layers.x.runs(tmp.csr$MAP$admix.proportions,csr$MAP$admix.proportions,csr1.order)
	}
	if(k > 2){
		tmp.csr <- output.list.nsp[[k-1]][[1]]
		csr1.order <- match.layers.x.runs(tmp.csr$MAP$admix.proportions,csr$MAP$admix.proportions,csr1.order)
	}
	if(is.null(csr1.order)){
		csr1.order <- 1:k
	}
	laycon.nsp[1:k,k] <- calculate.layer.contribution(csr,data.block,csr1.order)
}

pdf(file="sims/simK1_laycon_barplots.pdf",width=8,height=4,pointsize=14)
	#quartz(width=8,height=4,pointsize=14)
	par(mfrow=c(1,2),mar=c(4,4,3,0.5))
	barplot(laycon.sp,	
			col=layer.colors,
			xlab="",ylab="layer contribution")
	barplot(laycon.nsp,
			col=layer.colors,
			xlab="",ylab="")
			mtext(side=1,text="number of layers",padj=4.25,adj=-1)
			mtext(side=3,text=bquote(paste("Layer contributions (true ",italic("K"),"=1)")),padj=-1.3,adj=12,font=2,cex=1.2)
dev.off()


################################
#K2
################################
pdf(file="sims/simK2_std_xval.pdf",width=10,height=5,pointsize=14)
#	quartz(width=10,height=5)
	par(mfrow=c(1,2),mar=c(4,5,4,2))
	plot.sim.xvals(dir="../../sims/analyses/conStruct/cross_validation/K_2",n.reps=10,K=7,simK=2,y.lim=c(-10,0))
dev.off()


load("../../sims/sim_data/K_2/sim.dataset.Robj")
freq.data <- conStruct:::process.freq.data(sim.dataset$freq.data$freqs)
data.block <- conStruct:::make.data.block(K = 2,
										  freq.data = freq.data,
										  coords = sim.dataset$coords,
										  spatial = TRUE,
										  geoDist = fields::rdist(sim.dataset$coords))

output.list.sp <- vector("list",7)
for(k in 1:7){
	load(sprintf("../../sims/analyses/conStruct/runs/K_2/simK2_K%s_sp_conStruct.results.Robj",k))
	conStruct.results <- cluster.2.layer(conStruct.results)
	for(j in 1:k){
		names(conStruct.results[[1]]$MAP$layer.params[[j]])[4] <- "phi"
	}
	output.list.sp[[k]] <- conStruct.results
}

output.list.nsp <- vector("list",7)
for(k in 1:7){
	load(sprintf("../../sims/analyses/conStruct/runs/K_2/simK2_K%s_nsp_conStruct.results.Robj",k))
	conStruct.results <- cluster.2.layer(conStruct.results)
	for(j in 1:k){
		names(conStruct.results[[1]]$MAP$layer.params[[j]])[4] <- "phi"
	}
	output.list.nsp[[k]] <- conStruct.results
}

pdf(file="sims/simK2_nsp_pies.pdf",width=8,height=6.3,pointsize=14)
	K <- 7
	output.list <- output.list.nsp
	radii <- 1.7
	trueK <- 2
	par(mfrow=c(2,3),oma=c(1,0,3.5,0))
	layer.colors <- c("blue","red","goldenrod1","forestgreen","darkorchid1","deepskyblue","darkorange1","seagreen2","yellow1","black")
	csr1.order <- NULL
	for(k in 2:5){
		data.block$K <- k
		csr <- output.list[[k]][[1]]
		if(k > 2){
			tmp.csr <- output.list[[k-1]][[1]]
			csr1.order <- conStruct:::match.layers.x.runs(tmp.csr$MAP$admix.proportions,csr$MAP$admix.proportions,csr1.order)
			if(k ==4){
				tmp.4.csr.order <- csr1.order
			}
		}
		if(is.null(csr1.order)){
			csr1.order <- 1:k
		}
		make.admix.pie.plot(csr$MAP$admix.proportions[,csr1.order],data.block$coords,layer.colors=layer.colors,radii=radii,x.lim=c(2.5,8.5),y.lim=c(2.5,8.5),mar=c(4,2,1,2))
		if(k==3){
			mtext(side=3,text=bquote(paste("True ",italic("K")," = ",.(trueK))),cex=1.5,padj=-0.7)
		}
			mtext(side=1,text=bquote(paste("(",.(letters[k-1]),") ",italic("K")," = ",.(k))),padj=2.7,adj=0.4)
	}
	k <- 6
	data.block$K <- k
		csr <- output.list[[k]][[1]]
		tmp.csr <- output.list[[k-2]][[1]]
		csr1.order <- conStruct:::match.layers.x.runs(tmp.csr$MAP$admix.proportions,csr$MAP$admix.proportions,tmp.4.csr.order)
		make.admix.pie.plot(csr$MAP$admix.proportions[,csr1.order],data.block$coords,layer.colors=layer.colors,radii=radii,x.lim=c(2.5,8.5),y.lim=c(2.5,8.5),mar=c(4,2,1,2))
		mtext(side=1,text=bquote(paste("(",.(letters[k-1]),") ",italic("K")," = ",.(k))),padj=2.7,adj=0.4)
	k <- 7
	data.block$K <- k
		csr <- output.list[[k]][[1]]
		tmp.csr <- output.list[[k-1]][[1]]
		csr1.order <- conStruct:::match.layers.x.runs(tmp.csr$MAP$admix.proportions,csr$MAP$admix.proportions,csr1.order)
		make.admix.pie.plot(csr$MAP$admix.proportions[,csr1.order],data.block$coords,layer.colors=layer.colors,radii=radii,x.lim=c(2.5,8.5),y.lim=c(2.5,8.5),mar=c(4,2,1,2))
		mtext(side=1,text=bquote(paste("(",.(letters[k-1]),") ",italic("K")," = ",.(k))),padj=2.7,adj=0.4)
dev.off()

pdf(file="sims/simK2_sp_pies.pdf",width=8,height=6.3,pointsize=14)
	par(mfrow=c(2,3),oma=c(1,0,3.5,0))
	plot.sim.pies.multipanel(data.block = data.block,
							 K = 7,
							 output.list = output.list.sp,
							 radii = 1.7,
							 mar = c(4,2,1,2),trueK=2)
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
		csr1.order <- match.layers.x.runs(tmp.csr$MAP$admix.proportions,csr$MAP$admix.proportions,csr1.order)
	}
	if(is.null(csr1.order)){
		csr1.order <- 1:k
	}
	laycon.sp[1:k,k] <- calculate.layer.contribution(csr,data.block,csr1.order)
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
		csr1.order <- match.layers.x.runs(tmp.csr$MAP$admix.proportions,csr$MAP$admix.proportions,csr1.order)
	}
	if(k > 2){
		tmp.csr <- output.list.nsp[[k-1]][[1]]
		csr1.order <- match.layers.x.runs(tmp.csr$MAP$admix.proportions,csr$MAP$admix.proportions,csr1.order)
	}
	if(is.null(csr1.order)){
		csr1.order <- 1:k
	}
	laycon.nsp[1:k,k] <- calculate.layer.contribution(csr,data.block,csr1.order)
}

pdf(file="sims/simK2_laycon_barplots.pdf",width=8,height=4,pointsize=14)
	par(mfrow=c(1,2),mar=c(4,4,3,0.5))
	barplot(laycon.sp,	
			col=layer.colors,
			xlab="",ylab="layer contribution")
	barplot(laycon.nsp,
			col=layer.colors,
			xlab="",ylab="")
			mtext(side=1,text="number of layers",padj=4.25,adj=-1)
			mtext(side=3,text=bquote(paste("Layer contributions (true ",italic("K"),"=2)")),padj=-1.3,adj=12,font=2,cex=1.2)
dev.off()

pdf(file="sims/simK2_adprop_fit.pdf",width=6,height=6,pointsize=14)
	viz.admix.results(sim.admix.props = sim.dataset$admix.list$w,
				 	 conStruct.results = output.list.sp[[2]][[1]],
				 	 layer.order=c(2,1))
dev.off()

################################
#K3
################################
pdf(file="sims/simK3_std_xval.pdf",width=10,height=5,pointsize=14)
#	quartz(width=10,height=5)
	par(mfrow=c(1,2),mar=c(4,5,4,2))
	plot.sim.xvals(dir="../../sims/analyses/conStruct/cross_validation/K_3",n.reps=10,K=7,simK=3,y.lim=c(-10,0))
dev.off()


load("../../sims/sim_data/K_3/sim.dataset.Robj")
freq.data <- conStruct:::process.freq.data(sim.dataset$freq.data$freqs)
data.block <- conStruct:::make.data.block(K = 1,
										  freq.data = freq.data,
										  coords = sim.dataset$coords,
										  spatial = TRUE,
										  geoDist = fields::rdist(sim.dataset$coords))

output.list.sp <- vector("list",7)
for(k in 1:7){
	load(sprintf("../../sims/analyses/conStruct/runs/K_3/simK3_K%s_sp_conStruct.results.Robj",k))
	conStruct.results <- cluster.2.layer(conStruct.results)
	for(j in 1:k){
		names(conStruct.results[[1]]$MAP$layer.params[[j]])[4] <- "phi"
	}
	output.list.sp[[k]] <- conStruct.results
}

output.list.nsp <- vector("list",7)
for(k in 1:7){
	load(sprintf("../../sims/analyses/conStruct/runs/K_3/simK3_K%s_nsp_conStruct.results.Robj",k))
	conStruct.results <- cluster.2.layer(conStruct.results)
	for(j in 1:k){
		names(conStruct.results[[1]]$MAP$layer.params[[j]])[4] <- "phi"
	}
	output.list.nsp[[k]] <- conStruct.results
}


pdf(file="sims/simK3_nsp_pies.pdf",width=8,height=6.3,pointsize=14)
	par(mfrow=c(2,3),oma=c(1,0,3.5,0))
	plot.sim.pies.multipanel(data.block = data.block,
							 K = 7,
							 output.list = output.list.nsp,
							 radii = 1.7,
							 mar = c(4,2,1,2),trueK=3)
dev.off()

pdf(file="sims/simK3_sp_pies.pdf",width=8,height=6.3,pointsize=14)
	par(mfrow=c(2,3),oma=c(1,0,3.5,0))
	plot.sim.pies.multipanel(data.block = data.block,
							 K = 7,
							 output.list = output.list.sp,
							 radii = 1.7,
							 mar = c(4,2,1,2),trueK=3)
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
		csr1.order <- match.layers.x.runs(tmp.csr$MAP$admix.proportions,csr$MAP$admix.proportions,csr1.order)
	}
	if(is.null(csr1.order)){
		csr1.order <- 1:k
	}
	laycon.sp[1:k,k] <- calculate.layer.contribution(csr,data.block,csr1.order)
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
		csr1.order <- match.layers.x.runs(tmp.csr$MAP$admix.proportions,csr$MAP$admix.proportions,csr1.order)
	}
	if(k > 2){
		tmp.csr <- output.list.nsp[[k-1]][[1]]
		csr1.order <- match.layers.x.runs(tmp.csr$MAP$admix.proportions,csr$MAP$admix.proportions,csr1.order)
	}
	if(is.null(csr1.order)){
		csr1.order <- 1:k
	}
	laycon.nsp[1:k,k] <- calculate.layer.contribution(csr,data.block,csr1.order)
}

pdf(file="sims/simK3_laycon_barplots.pdf",width=8,height=4,pointsize=14)
	par(mfrow=c(1,2),mar=c(4,4,3,0.5))
	barplot(laycon.sp,	
			col=layer.colors,
			xlab="",ylab="layer contribution")
	barplot(laycon.nsp,
			col=layer.colors,
			xlab="",ylab="")
			mtext(side=1,text="number of layers",padj=4.25,adj=-1)
			mtext(side=3,text=bquote(paste("Layer contributions (true ",italic("K"),"=3)")),padj=-1.3,adj=12,font=2,cex=1.2)
dev.off()

pdf(file="sims/simK3_adprop_fit.pdf",width=6,height=6,pointsize=14)
viz.admix.results(sim.admix.props = sim.dataset$admix.list$w,
				  conStruct.results = output.list.sp[[3]][[1]],
				  layer.order=c(3,2,1))
dev.off()


################################
# K1 vs. K2 vs. K3
################################

CV.error1 <- get.CV.error(Rout.file="../../sims/analyses/admixture/datasets/simK1/exe.admixture.Rout")
CV.error2 <- get.CV.error(Rout.file="../../sims/analyses/admixture/datasets/simK2/exe.admixture.Rout")
CV.error3 <- get.CV.error(Rout.file="../../sims/analyses/admixture/datasets/simK3/exe.admixture.Rout")


pdf(file="sims/sim_xvals.pdf",width=14,height=5,pointsize=20)
#	quartz(width=14,height=5)
#K1
n.reps <- 10
K <- 7
for(n in 1:n.reps){
	load(sprintf("../../sims/analyses/conStruct/cross_validation/K_1/simK1_rep%s_test.lnl.Robj",n))
	assign(paste0("tl",n),test.lnl)
}
x.vals <- lapply(1:n.reps,function(n){get(sprintf("tl%s",n))})
x.vals.std <- lapply(x.vals,conStruct:::standardize.xvals)
xval.CIs <- conStruct:::get.xval.CIs(x.vals.std,K)
	par(mfrow=c(1,3),mar=c(4.5,4.5,4,1))
	plot.xval.CIs(xval.CIs,K,xlim=c(0.8,7.2),jitter=0.15)
		mtext("Predictive accuracy",side=2,padj=-3)
	legend(x="bottomright",pch=c(19,19,5),
		col=c("blue","green","orangered"),
			legend=c("spatial","nonspatial",expression("ADMIXTURE best "*italic("K"))),cex=1)
	mtext(bquote(paste("True ",italic("K")," = 1")),side=3,adj=0.5,padj=-1.5,font=2,cex=1.2)
		points(which.min(CV.error1)+0.15,xval.CIs$nsp.means[which.min(CV.error1)],pch=5,col="orangered",cex=1.5)
#K2
n.reps <- 10
K <- 7
for(n in 1:n.reps){
	load(sprintf("../../sims/analyses/conStruct/cross_validation/K_2/simK2_rep%s_test.lnl.Robj",n))
	assign(paste0("tl",n),test.lnl)
}
x.vals <- lapply(1:n.reps,function(n){get(sprintf("tl%s",n))})
x.vals.std <- lapply(x.vals,conStruct:::standardize.xvals)
xval.CIs <- conStruct:::get.xval.CIs(x.vals.std,K)
	plot.xval.CIs(xval.CIs,K,xlim=c(0.8,7.2),jitter=0.15)
		points(which.min(CV.error2)+0.15,xval.CIs$nsp.means[which.min(CV.error2)],pch=5,col="orangered",cex=1.5)
	rect(0.7,-700,7.3,400,lty=2)
	arrows(0.7,-700,2.5,-3e3,lty=1,length=0.15)
	arrows(7.2,-700,6.6,-3e3,lty=1,length=0.15)
	TeachingDemos::subplot(fun = {
						plot.xval.CIs(xval.CIs,K,k.range=c(1:7),xlim=c(0.7,7.2),ylim=c(-275,0),axes=FALSE,cex=1,jitter=0.1)
							axis(1,at=1:7,labels=c(1,"","","","","",7),cex.axis=0.8,lty=1)
							axis(2,at=seq(-275,0,length.out=6),labels=c(-275,"","","","",0),cex.axis=0.8,lty=1)
							box(lwd=1.2,lty=2)
							box(lwd=1.2,lty=1,bty="l")
							points(which.min(CV.error2)+0.15,xval.CIs$nsp.means[which.min(CV.error2)],pch=5,col="orangered",cex=1)
						},
					x=c(2.5,6.6),y=c(-12.5e3,-3e3))
	mtext(bquote(paste("True ",italic("K")," = 2")),side=3,adj=0.5,padj=-1.5,font=2,cex=1.2)
	mtext("number of layers",side=1,padj=3.4)
#K3
n.reps <- 10
K <- 7
for(n in 1:n.reps){
	load(sprintf("../../sims/analyses/conStruct/cross_validation/K_3/simK3_rep%s_test.lnl.Robj",n))
	assign(paste0("tl",n),test.lnl)
}
x.vals <- lapply(1:n.reps,function(n){get(sprintf("tl%s",n))})
x.vals.std <- lapply(x.vals,conStruct:::standardize.xvals)
xval.CIs <- conStruct:::get.xval.CIs(x.vals.std,K)
	plot.xval.CIs(xval.CIs,K,xlim=c(0.8,7.2),jitter=0.15)
		points(which.min(CV.error3)+0.15,xval.CIs$nsp.means[which.min(CV.error3)],pch=5,col="orangered",cex=1.5)
	rect(0.7,-680,7.3,470,lty=2)
	arrows(0.7,-680,3.2,-4e3,lty=1,length=0.1)
	arrows(7.2,-700,6.6,-4e3,lty=1,length=0.1)
	TeachingDemos::subplot(fun = {
						plot.xval.CIs(xval.CIs,K,k.range=c(1:7),ylim=c(-140,0),xlim=c(0.7,7.3),axes=FALSE,cex=1,jitter=0.1)
							axis(1,at=1:7,labels=c(1,"","","","","",7),cex.axis=0.8,lty=1)
							axis(2,at=seq(-140,0,length.out=6),labels=c(-140,"","","","",0),cex.axis=0.8,lty=1)
							box(lwd=1.2,lty=2)
							box(lwd=1.2,lty=1,bty="l")
							points(which.min(CV.error3)+0.15,xval.CIs$nsp.means[which.min(CV.error3)],pch=5,col="orangered",cex=1)
						},
					x=c(3.2,6.6),y=c(-1.5e4,-4e3))
	mtext(bquote(paste("True ",italic("K")," = 3")),side=3,adj=0.5,padj=-1.5,font=2,cex=1.2)
dev.off()