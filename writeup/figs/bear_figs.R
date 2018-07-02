source("fig_funcs.R")

n.reps <- 10
K <- 7
for(n in 1:n.reps){
	load(sprintf("../../data/bears/conStruct/xvalidation/bear_rep%s_test.lnl.Robj",n,n))
	assign(paste0("tl",n),test.lnl)
}
x.vals <- lapply(1:n.reps,function(n){get(sprintf("tl%s",n))})
x.vals.std <- lapply(x.vals,conStruct:::standardize.xvals)
xval.CIs <- conStruct:::get.xval.CIs(x.vals.std,K)
pdf(file="bears/bear_std_xval.pdf",width=12,height=4.25,pointsize=14)
#	quartz(width=12,height=5)
	par(mfrow=c(1,3),mar=c(3.75,5,4,1))
	plot.xval.CIs(xval.CIs,K,jitter=0.1,xlim=c(0.75,7.25))
		legend(x="bottomright",pch=19,col=c("blue","green"),legend=c("spatial","nonspatial"))
	mtext("Predictive accuracy",side=2,padj=-3,font=2)
	mtext("conStruct",side=3,font=2)
	plot.xval.CIs(xval.CIs,K,ylim=c(-200,0),jitter=0.1,xlim=c(3.75,7.25),xaxt='n')
		axis(1,at=c(4:7))
		legend(x="bottomright",pch=c(19,NA),lty=c(NA,1),lwd=c(NA,2),col=c(1,adjustcolor(1,0.8)),legend=c("mean","95% CI"))
	mtext("number of layers",side=1,font=2,padj=2.75)
	mtext("Cross-validation results (Bears)",side=3,padj=-2.2,font=2,cex=1.2)
	mtext("non-spatial conStruct",side=3,font=2)
	mtext("Predictive accuracy",side=2,padj=-3,font=2)
	CV.error <- get.CV.error("../../data/bears/admixture/exe.admixture.Rout")
	plot(CV.error,pch=19,col="orangered",
			xlab="",ylab="",cex=2,ylim=rev(range(CV.error)))
	best <- which.min(CV.error)
	points(best,CV.error[best],pch=8,cex=2)
	legend(x="bottomright",pch=8,legend="model with minimum CV error")
	mtext("ADMIXTURE",side=3,font=2)
	mtext("Cross-Validation Error",side=2,padj=-3,font=2)
dev.off()

library(maps)
load("../../data/bears/data/bear.dataset.Robj")
pdf(file="bears/bear_sampling_map.pdf",width=6,height=4,pointsize=13)
	#quartz(width=6,height=4,pointsize=13)
	par(mar=c(1,1.5,0,0),oma=c(0,0,0,0))
	map(xlim = range(bear.dataset$sample.coords[,1]) + c(-5,5), ylim = range(bear.dataset$sample.coords[,2])+c(-2,2), col="gray",lforce="e")
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
	load(sprintf("../../data/bears/conStruct/runs/bearsK%s_sp_conStruct.results.Robj",k))
	conStruct.results <- cluster.2.layer(conStruct.results)
	for(j in 1:k){
		names(conStruct.results[[1]]$MAP$layer.params[[j]])[4] <- "phi"
	}
	output.list.sp[[k]] <- conStruct.results
}

output.list.nsp <- vector("list",7)
for(k in 1:7){
	load(sprintf("../../data/bears/conStruct/runs/bearsK%s_nsp_conStruct.results.Robj",k))
	conStruct.results <- cluster.2.layer(conStruct.results)
	for(j in 1:k){
		names(conStruct.results[[1]]$MAP$layer.params[[j]])[4] <- "phi"
	}
	output.list.nsp[[k]] <- conStruct.results
}

pdf(file="bears/bears_sp_vs_nsp.pdf",width=22.5,height=7.5,pointsize=14)
	layout(cbind(matrix(c(rep(1,10),rep(2,15)),nrow=5,ncol=5,byrow=TRUE),
				 matrix(2+c(rep(1,10),rep(2,15)),nrow=5,ncol=5,byrow=TRUE),
				 matrix(4+c(rep(1,10),rep(2,15)),nrow=5,ncol=5,byrow=TRUE)))
	par(mar=c(5,5,1,1))
	make.bear.redux.result.plot.multipanel1(admix.proportions = output.list.sp[[3]][[1]]$MAP$admix.proportions,
										coords = bear.dataset$sample.coords,
										lump.dist = 200,
										layer.colors = layer.colors[order(c(1,2,3))],
										layer.order=c(1,2,3))
		mtext(side=1,text=bquote(paste("(",.(letters[1]),") ",italic("K")," = ",.(3))),padj=-2.3,adj=0.03,cex=1.3)
		mtext(side=1,text="conStruct (spatial)",padj=-1.6,adj=0.12,cex=1.3)
#		mtext(side=1,text=bquote(paste("(",.(letters[1]),") ",italic("K")," = ",.(3),"  (spatial)")),padj=-1.5,adj=0.03,cex=1.3)
	par(xpd=FALSE)
	make.bear.redux.result.plot.multipanel1(admix.proportions = output.list.nsp[[3]][[1]]$MAP$admix.proportions[,c(3,1,2)],
										coords = bear.dataset$sample.coords,
										lump.dist = 200,
										layer.colors = layer.colors[order(c(1,2,3))],
										layer.order=c(1,2,3))
		mtext(side=1,text=bquote(paste("(",.(letters[2]),") ",italic("K")," = ",.(3))),padj=-2.3,adj=0.03,cex=1.3)
		mtext(side=1,text="conStruct (nonspatial)",padj=-1.6,adj=0.12,cex=1.3)
#		mtext(side=1,text=bquote(paste("(",.(letters[2]),") ",italic("K")," = ",.(3)," (nonspatial)")),padj=-1.5,adj=0.03,cex=1.3)
	par(xpd=FALSE)
		w <- as.matrix(read.table("../../data/bears/admixture/bears.3.Q",stringsAsFactors=FALSE))
		make.bear.redux.result.plot.multipanel1(admix.proportions = w[,c(3,2,1)],
												coords = bear.dataset$sample.coords,
												lump.dist = 200,
												layer.colors = layer.colors,
												layer.order=NULL)
		mtext(side=1,text=bquote(paste("(",.(letters[3]),") ",italic("K")," = ",.(3))),padj=-2.3,adj=0.03,cex=1.3)
		mtext(side=1,text="ADMIXTURE",padj=-1.6,adj=0.1,cex=1.3)
#		mtext(side=1,text=bquote(paste("(",.(letters[3]),") ",italic("K")," = ",.(3)," (ADMIXTURE)")),padj=-1.5,adj=0.03,cex=1.3)
dev.off()

pdf(file="bears/bear_sp_results.pdf",width=15,height=10)
	layout(
		rbind(
			Reduce("cbind",lapply(1:3,function(x){
				matrix(((x-1)*2)+c(rep(1,10),rep(2,15)),
					nrow=5,ncol=5,byrow=TRUE)})),
			Reduce("cbind",lapply(4:6,function(x){
				matrix(((x-1)*2)+c(rep(1,10),rep(2,15)),
					nrow=5,ncol=5,byrow=TRUE)}))
			))
	make.bear.redux.result.plot.multipanel2(output.list = output.list.sp,
											coords = data.block$coords,
											lump.dist = 200,
											layer.colors,csr1.order=c(2,1))
dev.off()

pdf(file="bears/bear_nsp_results.pdf",width=15,height=10)
	layout(
		rbind(
			Reduce("cbind",lapply(1:3,function(x){
				matrix(((x-1)*2)+c(rep(1,10),rep(2,15)),
					nrow=5,ncol=5,byrow=TRUE)})),
			Reduce("cbind",lapply(4:6,function(x){
				matrix(((x-1)*2)+c(rep(1,10),rep(2,15)),
					nrow=5,ncol=5,byrow=TRUE)}))
			))
	make.bear.redux.result.plot.multipanel2(output.list = output.list.nsp,
											coords = data.block$coords,
											lump.dist = 200,
											layer.colors,csr1.order=NULL)
dev.off()


pdf(file="bears/bear_sp_layer_covs.pdf",width=12,height=8,pointsize=14)
	#quartz(width=12,height=8)
	layout(matrix(c(1:6),nrow=2,ncol=3,byrow=TRUE))
	par(mar=c(4,5,3,2),oma=c(3,3,3,1))	
	plot.K.layer.curves(K=1:6,data.block=data.block,output.list= output.list.sp,col.mat1=NULL,col.mat2=NULL,y.range=c(0.15,0.25))
	mtext(text="geographic distance (mi)",side=1,font=2,cex.axis=2,padj=4.5,adj=-6.2)
	mtext(text="allele frequency covariance",side=2,font=2,cex.axis=2,padj=-59.5,adj=20)
dev.off()


pdf(file="bears/bear_nsp_layer_covs.pdf",width=12,height=8,pointsize=14)
	#quartz(width=12,height=8)
	layout(matrix(c(1:6),nrow=2,ncol=3,byrow=TRUE))
	par(mar=c(4,5,3,2),oma=c(3,3,3,1))	
	plot.K.layer.curves(K=1:6,data.block=data.block,output.list= output.list.nsp,col.mat1=NULL,col.mat2=NULL,output.list.sp=output.list.sp)
	mtext(text="geographic distance (mi)",side=1,font=2,cex.axis=2,padj=4.5,adj=-6.2)
	mtext(text="allele frequency covariance",side=2,font=2,cex.axis=2,padj=-59.5,adj=20)
dev.off()

pdf(file="bears/bear_sp_colorized_cov.pdf",width=8,height=8,pointsize=14)
	sd.D <- sd(fields::rdist.earth(data.block$coords))
	lon.cols <- terrain.colors(data.block$N)[as.numeric(cut(data.block$coords[,1],data.block$N))]
	col.mats <- make.col.mats(lon.cols)
	data.block$K <- 2
	plot(data.block$geoDist,data.block$obsCov,xlab="",ylab="",xaxt='n',
			pch=21,col=col.mats[[2]],bg=col.mats[[1]])
		axis(side=1,at=c(seq(min(data.block$geoDist),max(data.block$geoDist),length.out=5)),
				round(seq(min(data.block$geoDist),
							max(data.block$geoDist),length.out=5) * sd.D,0))
		mtext(text="geographic distance (mi)",side=1,font=2,cex.axis=2,padj=3.7,adj=0.5)
		mtext(text="allele frequency covariance",side=2,font=2,cex.axis=2,padj=-3.7,adj=0.5)
	plot.layer.curves(data.block=data.block,conStruct.results=output.list.sp[[2]][[1]],
						layer.cols=NULL,add=TRUE,col.mat1=NULL,col.mat2=NULL)
	TeachingDemos::subplot(fun = {
						plot(0,xlim=range(data.block$coords[,1]) + c(-5,5),
							  ylim = range(data.block$coords[,2])+c(-2,2), 
							  type='n',yaxt='n',xaxt='n',xlab="",ylab="")
							map(xlim = range(data.block$coords[,1]) + c(-5,5), 
								ylim = range(data.block$coords[,2])+c(-2,2), 
								col="gray",mar=c(1,1,1,1),add=TRUE)
							points(data.block$coords,bg=lon.cols,pch=21,col=1,cex=2)
							box(lwd=1.1)
						},
					x=c(1.3e3/sd.D,3.8e3/sd.D),y=c(0.205,0.245))
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

pdf(file="bears/bears_laycon_barplots.pdf",width=8,height=4,pointsize=14)
	par(mfrow=c(1,2),mar=c(4,4,3,0.5))
	barplot(laycon.sp,	
			col=layer.colors,
			xlab="",ylab="layer contribution")
	barplot(laycon.nsp,
			col=layer.colors,
			xlab="",ylab="")
			mtext(side=1,text="number of layers",padj=4.25,adj=-1)
			mtext(side=3,text="Layer contributions (Bears)",padj=-2.25,adj=18,font=2,cex=1.2)
dev.off()


pdf(file="bears/bear_admixture_results.pdf",width=15,height=10,pointsize=14)
	# layout(Reduce("cbind",lapply(1:3,function(x){
			# matrix(((x-1)*2)+c(rep(1,10),rep(2,15)),
					# nrow=5,ncol=5,byrow=TRUE)})))
	layout(
		rbind(
			Reduce("cbind",lapply(1:3,function(x){
				matrix(((x-1)*2)+c(rep(1,10),rep(2,15)),
					nrow=5,ncol=5,byrow=TRUE)})),
			Reduce("cbind",lapply(4:6,function(x){
				matrix(((x-1)*2)+c(rep(1,10),rep(2,15)),
					nrow=5,ncol=5,byrow=TRUE)}))
			))
	for(k in 2:7){
		w <- as.matrix(read.table(sprintf("../../data/bears/admixture/bears.%s.Q",k),stringsAsFactors=FALSE))
		if(k <= 2){
			csr1.order <- c(2,1)
		}
		if(k > 2){
			tmp.w <- as.matrix(read.table(sprintf("../../data/bears/admixture/bears.%s.Q",k-1),stringsAsFactors=FALSE))
			csr1.order <- match.layers.x.runs(tmp.w,w,csr1.order)
		}
		if(is.null(csr1.order)){
			csr1.order <- 1:k
		}
		make.bear.redux.result.plot.multipanel1(admix.proportions = w[,csr1.order],
												coords = bear.dataset$sample.coords,
												lump.dist = 200,
												layer.colors = layer.colors,
												layer.order=NULL)
		mtext(side=1,text=bquote(paste("(",.(letters[k-1]),") ",italic("K")," = ",.(k))),padj=0.8,adj=0.4)
	}
dev.off()

pdf(file="bears/bear_admixture_CVerror.pdf",width=7,height=7,pointsize=14)
	CV.error <- get.CV.error("../../data/bears/admixture/exe.admixture.Rout")
	plot(CV.error,pch=19,col="orangered",
			xlab="number of clusters",ylab="Cross-Validation Error",cex=2,
			main="Bear ADMIXTURE cross-validation results")
	best <- which.min(CV.error)
	points(best,CV.error[best],pch=8,cex=2)
	legend(x="topright",pch=8,legend="model with minimum CV error")
dev.off()