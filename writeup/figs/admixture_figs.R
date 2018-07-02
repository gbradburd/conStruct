get.CV.error <- function(Rout.file){
	log <- scan(Rout.file,what="character",sep="\n")
	CV.error <- as.numeric(
					unlist(
						lapply(
							strsplit(
								log[grepl("CV error",log)],
								": "),
						"[[",2)
					)
				)
	return(CV.error)
}

collapse.ind.Q <- function(Q,pop.vec){
	n.pops <- length(unique(pop.vec))
	admix.props <- matrix(NA,n.pops,ncol(Q))
	for(i in 1:n.pops){
		admix.props[i,] <- colMeans(Q[which(pop.vec==i),,drop=FALSE])
	}
	return(admix.props)
}

pie.plot <- function(Q.file,coords,pop.vec,layer.colors){
	Q <- read.table(Q.file,stringsAsFactors=FALSE)
	admix.props <- collapse.ind.Q(Q,pop.vec)
	K <- ncol(admix.props)
	N <- length(unique(pop.vec))
	layer.names <- paste0("layer_",1:K)
	sample.names <- paste0("sample_", 1:N)
	color.tab <- caroline::nv(c(layer.colors[1:K]), layer.names)
	pie.list <- lapply(1:N, function(i) {
	    caroline::nv(admix.props[i, ], layer.names)
	})
	names(pie.list) <- sample.names
	caroline::pies(pie.list, x0 = coords[, 1], 
	    y0 = coords[, 2], color.table = color.tab, 
	    border = "black", radii = 3, xlab = "", ylab = "", 
	    lty = 1, density = NULL,
	    xlim = range(coords[,1]) + c(-0.5,0.5),
	    ylim = range(coords[,2]) + c(-0.5,0.5))
	box(lwd = 2)
}

get.adprops <- function(Q.file,pop.vec){
	Q <- read.table(Q.file,stringsAsFactors=FALSE)
	admix.props <- collapse.ind.Q(Q,pop.vec)
	return(admix.props)
}

get.laycols <- function(K,q.file,N){
	layer.colors <- c("blue", "red", "goldenrod1", "forestgreen", 
    		          "darkorchid1", "deepskyblue", "darkorange1", "seagreen2", 
    		          "yellow1", "black")
	fastStr.admix.props <- lapply(1:K,
								function(k){
									get.adprops(Q.file = sprintf(q.file,k),
												pop.vec = unlist(lapply(1:N,function(n){rep(n,10)})))
								})
	layord <- NULL
	laycols <- list()
	for(k in 2:K){
		if(k > 2){
			layord <- match.layers.x.runs(fastStr.admix.props[[k-1]],
									  	  fastStr.admix.props[[k]],
									  	  layord)
		}
		if(is.null(layord)){
			layord <- 1:k
		}
		laycols[[k-1]] <- layer.colors[order(layord)]
	}
	return(laycols)
}

load("../../sims/sim_data/K_1/sim.dataset.Robj")
laycols <- get.laycols(K=7,q.file="../../sims/analyses/admixture/datasets/simK1/simK1.%s.Q",N=sim.dataset$N)


pdf(file="admixture/admixture_simK1_pies.pdf",width=12,height=7.75,pointsize=13)
	par(mfrow=c(2,3),mar=c(4.5,4,2.5,4))
	for(k in 2:7){
		pie.plot(Q.file = sprintf("../../sims/analyses/admixture/datasets/simK1/simK1.%s.Q",k),
			 coords = sim.dataset$coords,
			 pop.vec = unlist(lapply(1:sim.dataset$N,function(n){rep(n,10)})),
			 layer.colors = laycols[[k-1]])
		mtext(side=3,font=2,text=bquote(paste(italic("K"),"=",.(k))))
		if(k == 6){
			mtext(side=1,font=2,text=bquote(paste("true ",italic("K")," = ",1)),cex=1.2,padj=2.7)
		}
	}
dev.off()

ad.Qs <- lapply(2:7,
			function(k){
				pop.vec <- unlist(lapply(1:sim.dataset$N,function(n){rep(n,10)}))
				Q.file <- sprintf("../../sims/analyses/admixture/datasets/simK1/simK1.%s.Q",k)
				Q <- read.table(Q.file,stringsAsFactors=FALSE)
				admix.props <- collapse.ind.Q(Q,pop.vec)
				return(admix.props)
			})

for(k in 1:6){
	pdf(file=sprintf("admixture/simK1_%s_struct.pdf",k+1),width=10,height=4.5)
		make.structure.plot(ad.Qs[[k]],layer.colors=laycols[[k]])
	dev.off()
	pdf(file=sprintf("admixture/simK1_%s_pies.pdf",k+1),width=6,height=6)
		make.admix.pie.plot(ad.Qs[[k]],sim.dataset$coords,layer.colors=laycols[[k]],radii=6,x.lim=c(2.5,8.5),y.lim=c(2.5,8.5))
		mtext(side=3,font=2,text=bquote(paste(italic("K"),"=",.(k+1))),cex=2)
	dev.off()
}

load("../../sims/sim_data/K_2/sim.dataset.Robj")
laycols <- get.laycols(K=7,q.file="../../sims/analyses/admixture/datasets/simK2/simK2.%s.Q",N=sim.dataset$N)

pdf(file="admixture/admixture_simK2_pies.pdf",width=12,height=7.75,pointsize=13)
	par(mfrow=c(2,3),mar=c(4.5,4,2.5,4))
	for(k in 2:7){
		pie.plot(Q.file = sprintf("../../sims/analyses/admixture/datasets/simK2/simK2.%s.Q",k),
			 coords = sim.dataset$coords,
			 pop.vec = unlist(lapply(1:sim.dataset$N,function(n){rep(n,10)})),
			 layer.colors = laycols[[k-1]])
		mtext(side=3,font=2,text=bquote(paste(italic("K"),"=",.(k))))
		if(k == 6){
			mtext(side=1,font=2,text=bquote(paste("true ",italic("K")," = ",2)),cex=1.2,padj=2.7)
		}
	}
dev.off()
			
load("../../sims/sim_data/K_3/sim.dataset.Robj")
laycols <- get.laycols(K=7,q.file="../../sims/analyses/admixture/datasets/simK3/simK3.%s.Q",N=sim.dataset$N)

pdf(file="admixture/admixture_simK3_pies.pdf",width=12,height=7.75,pointsize=13)
	par(mfrow=c(2,3),mar=c(4.5,4,2.5,4))
	for(k in 2:7){
		pie.plot(Q.file = sprintf("../../sims/analyses/admixture/datasets/simK3/simK3.%s.Q",k),
			 coords = sim.dataset$coords,
			 pop.vec = unlist(lapply(1:sim.dataset$N,function(n){rep(n,10)})),
			 layer.colors = laycols[[k-1]])
		mtext(side=3,font=2,text=bquote(paste(italic("K"),"=",.(k))))
		if(k == 6){
			mtext(side=1,font=2,text=bquote(paste("true ",italic("K")," = ",3)),cex=1.2,padj=2.7)
		}
	}
dev.off()

CV.error1 <- get.CV.error(Rout.file="../../sims/analyses/admixture/datasets/simK1/exe.admixture.Rout")
CV.error2 <- get.CV.error(Rout.file="../../sims/analyses/admixture/datasets/simK2/exe.admixture.Rout")
CV.error3 <- get.CV.error(Rout.file="../../sims/analyses/admixture/datasets/simK3/exe.admixture.Rout")

pdf(file="admixture/sims_CVerror.pdf",width=15,height=5.5,pointsize=14)
	#quartz(width=15,height=5)
	par(mfrow=c(1,3),mar=c(4.5,4.5,4.9,2))
	plot(CV.error1,pch=19,col="orangered",
			xlab="",ylab="Cross-Validation Error",cex=2,cex.lab=1.7,cex.axis=1.5)
		best <- which.min(CV.error1)
		points(best,CV.error1[best],pch=8,cex=2)
	legend(x="topright",pch=8,legend="model with minimum CV error",cex=1.4)
	mtext(side=3,font=2,text=bquote(paste("True ",italic("K"),"= 1")),padj=-0.2,cex=1.1,adj=0)
	plot(CV.error2,pch=19,col="orangered",
			xlab="number of clusters",ylab="",cex=2,cex.lab=1.7,cex.axis=1.5)
		best <- which.min(CV.error2)
		points(best,CV.error2[best],pch=8,cex=2)
	mtext("ADMIXTURE cross-validation results for simulated data",cex=1.6,padj=-2.2)
	xl <- 1.9
	yb <- min(CV.error2)-0.01
	xr <- 7.1
	yt <- CV.error2[2]+0.01
	rect(xl,yb,xr,yt,lty=2)
		segments(xl,yt,2.2,0.3)
		segments(xr,yt,6.9,0.3)
		TeachingDemos::subplot(fun = {
							plot(CV.error2,pch=19,col="orangered",xlab="",ylab="",
									xlim=c(2,7),ylim=range(CV.error2[2:7]),cex=1.5,cex.axis=1.5);
							points(best,CV.error2[best],pch=8,cex=1.5)
						},
						x=c(2.2,6.9),y=c(0.3,0.65))
	mtext(side=3,font=2,text=bquote(paste("True ",italic("K"),"= 2")),padj=-0.2,cex=1.1,adj=0)
	plot(CV.error3,pch=19,col="orangered",
			xlab="",ylab="",cex=2,cex.lab=1.7,cex.axis=1.5)
		best <- which.min(CV.error3)
		points(best,CV.error3[best],pch=8,cex=2)
	xl <- 2.9
	yb <- min(CV.error3)-0.01
	xr <- 7.1
	yt <- CV.error3[3]+0.01
	rect(xl,yb,xr,yt,lty=2)
		segments(xl,yt,3.1,0.25)
		segments(xr,yt,6.9,0.25)
		TeachingDemos::subplot(fun = {
							plot(CV.error3,pch=19,col="orangered",xlab="",ylab="",
									xlim=c(3,7),ylim=range(CV.error3[3:7]),cex=1.5,cex.axis=1.5);
							points(best,CV.error3[best],pch=8,cex=1.5)
						},
						x=c(3.1,6.9),y=c(0.25,0.6))
	mtext(side=3,font=2,text=bquote(paste("True ",italic("K"),"= 3")),padj=-0.2,cex=1.1,adj=0)
dev.off()