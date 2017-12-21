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


load("~/Dropbox/conStruct/sims/cross_validation/K_1/sim.dataset.Robj")
laycols <- get.laycols(K=5,q.file="~/Dropbox/conStruct/sims/structure/datasets/simK1/simK1_estK.%s.meanQ",N=sim.dataset$N)


pdf(file="~/Dropbox/conStruct/writeup/figs/fastStr/fastStr_simK1_pies.pdf",width=12,height=4.5,pointsize=13)
	par(mfrow=c(2,3))
		pie.plot(Q.file = "~/Dropbox/conStruct/sims/structure/datasets/simK1/simK1_estK.2.meanQ",
				 coords = sim.dataset$coords,
				 pop.vec = unlist(lapply(1:sim.dataset$N,function(n){rep(n,10)})),
				 layer.colors = laycols[[1]])
			mtext(side=3,font=2,text=bquote(paste(italic("K"),"=",2)))
		pie.plot(Q.file = "~/Dropbox/conStruct/sims/structure/datasets/simK1/simK1_estK.3.meanQ",
				 coords = sim.dataset$coords,
				 pop.vec = unlist(lapply(1:sim.dataset$N,function(n){rep(n,10)})),
				 layer.colors = laycols[[2]])
			mtext(side=3,font=2,text=bquote(paste(italic("K"),"=",3)))
		pie.plot(Q.file = "~/Dropbox/conStruct/sims/structure/datasets/simK1/simK1_estK.4.meanQ",
				 coords = sim.dataset$coords,
				 pop.vec = unlist(lapply(1:sim.dataset$N,function(n){rep(n,10)})),
				 layer.colors = laycols[[3]])
			mtext(side=3,font=2,text=bquote(paste(italic("K"),"=",4)))
		pie.plot(Q.file = "~/Dropbox/conStruct/sims/structure/datasets/simK1/simK1_estK.5.meanQ",
				 coords = sim.dataset$coords,
				 pop.vec = unlist(lapply(1:sim.dataset$N,function(n){rep(n,10)})),
				 layer.colors = laycols[[4]])
			mtext(side=3,font=2,text=bquote(paste(italic("K"),"=",5)))
		# pie.plot(Q.file = "~/Dropbox/conStruct/sims/structure/datasets/simK1/simK1_estK.6.meanQ",
				 # coords = sim.dataset$coords,
				 # pop.vec = unlist(lapply(1:sim.dataset$N,function(n){rep(n,10)})),
				 # layer.colors = layer.colors[1:6])
			# mtext(side=3,font=2,text=bquote(paste(italic("K"),"=",6)))
			# mtext(side=1,font=2,text=bquote(paste("true ",italic("K")," = ",1)),cex=1.2,padj=2.7)
		# pie.plot(Q.file = "~/Dropbox/conStruct/sims/structure/datasets/simK1/simK1_estK.7.meanQ",
				 # coords = sim.dataset$coords,
				 # pop.vec = unlist(lapply(1:sim.dataset$N,function(n){rep(n,10)})),
				 # layer.colors = layer.colors[1:7])
			# mtext(side=3,font=2,text=bquote(paste(italic("K"),"=",7)))			
dev.off()


load("~/Dropbox/conStruct/sims/cross_validation/K_2/sim.dataset.Robj")
laycols <- get.laycols(K=6,q.file="~/Dropbox/conStruct/sims/structure/datasets/simK2/simK2_estK.%s.meanQ",N=sim.dataset$N)

pdf(file="~/Dropbox/conStruct/writeup/figs/fastStr/fastStr_simK2_pies.pdf",width=12,height=4.5,pointsize=13)
	par(mfrow=c(2,3))
		pie.plot(Q.file = "~/Dropbox/conStruct/sims/structure/datasets/simK2/simK2_estK.2.meanQ",
				 coords = sim.dataset$coords,
				 pop.vec = unlist(lapply(1:sim.dataset$N,function(n){rep(n,10)})),
				 layer.colors = laycols[[1]])
			mtext(side=3,font=2,text=bquote(paste(italic("K"),"=",2)))
		pie.plot(Q.file = "~/Dropbox/conStruct/sims/structure/datasets/simK2/simK2_estK.3.meanQ",
				 coords = sim.dataset$coords,
				 pop.vec = unlist(lapply(1:sim.dataset$N,function(n){rep(n,10)})),
				 layer.colors = laycols[[2]])
			mtext(side=3,font=2,text=bquote(paste(italic("K"),"=",3)))
		pie.plot(Q.file = "~/Dropbox/conStruct/sims/structure/datasets/simK2/simK2_estK.4.meanQ",
				 coords = sim.dataset$coords,
				 pop.vec = unlist(lapply(1:sim.dataset$N,function(n){rep(n,10)})),
				 layer.colors = laycols[[3]])
			mtext(side=3,font=2,text=bquote(paste(italic("K"),"=",4)))
		pie.plot(Q.file = "~/Dropbox/conStruct/sims/structure/datasets/simK2/simK2_estK.5.meanQ",
				 coords = sim.dataset$coords,
				 pop.vec = unlist(lapply(1:sim.dataset$N,function(n){rep(n,10)})),
				 layer.colors = laycols[[4]])
			mtext(side=3,font=2,text=bquote(paste(italic("K"),"=",5)))
		pie.plot(Q.file = "~/Dropbox/conStruct/sims/structure/datasets/simK2/simK2_estK.6.meanQ",
				 coords = sim.dataset$coords,
				 pop.vec = unlist(lapply(1:sim.dataset$N,function(n){rep(n,10)})),
				 layer.colors = laycols[[5]])
			mtext(side=3,font=2,text=bquote(paste(italic("K"),"=",6)))
			mtext(side=1,font=2,text=bquote(paste("true ",italic("K")," = ",2)),cex=1.2,padj=2.7)
		# pie.plot(Q.file = "~/Dropbox/conStruct/sims/structure/datasets/simK2/simK2_estK.7.meanQ",
				 # coords = sim.dataset$coords,
				 # pop.vec = unlist(lapply(1:sim.dataset$N,function(n){rep(n,10)})),
				 # layer.colors = layer.colors[1:7])
			# mtext(side=3,font=2,text=bquote(paste(italic("K"),"=",7)))
dev.off()
			
load("~/Dropbox/conStruct/sims/cross_validation/K_3/sim.dataset.Robj")
laycols <- get.laycols(K=7,q.file="~/Dropbox/conStruct/sims/structure/datasets/simK3/simK3_estK.%s.meanQ",N=sim.dataset$N)

pdf(file="~/Dropbox/conStruct/writeup/figs/fastStr/fastStr_simK3_pies.pdf",width=12,height=4.5,pointsize=13)
	par(mfrow=c(2,3))
		pie.plot(Q.file = "~/Dropbox/conStruct/sims/structure/datasets/simK3/simK3_estK.2.meanQ",
				 coords = sim.dataset$coords,
				 pop.vec = unlist(lapply(1:sim.dataset$N,function(n){rep(n,10)})),
				 layer.colors = laycols[[1]])
			mtext(side=3,font=2,text=bquote(paste(italic("K"),"=",2)))
		pie.plot(Q.file = "~/Dropbox/conStruct/sims/structure/datasets/simK3/simK3_estK.3.meanQ",
				 coords = sim.dataset$coords,
				 pop.vec = unlist(lapply(1:sim.dataset$N,function(n){rep(n,10)})),
				 layer.colors = laycols[[2]])
			mtext(side=3,font=2,text=bquote(paste(italic("K"),"=",3)))
		pie.plot(Q.file = "~/Dropbox/conStruct/sims/structure/datasets/simK3/simK3_estK.4.meanQ",
				 coords = sim.dataset$coords,
				 pop.vec = unlist(lapply(1:sim.dataset$N,function(n){rep(n,10)})),
				 layer.colors = laycols[[3]])
			mtext(side=3,font=2,text=bquote(paste(italic("K"),"=",4)))
		pie.plot(Q.file = "~/Dropbox/conStruct/sims/structure/datasets/simK3/simK3_estK.5.meanQ",
				 coords = sim.dataset$coords,
				 pop.vec = unlist(lapply(1:sim.dataset$N,function(n){rep(n,10)})),
				 layer.colors = laycols[[4]])
			mtext(side=3,font=2,text=bquote(paste(italic("K"),"=",5)))
		pie.plot(Q.file = "~/Dropbox/conStruct/sims/structure/datasets/simK3/simK3_estK.6.meanQ",
				 coords = sim.dataset$coords,
				 pop.vec = unlist(lapply(1:sim.dataset$N,function(n){rep(n,10)})),
				 layer.colors = laycols[[5]])
			mtext(side=3,font=2,text=bquote(paste(italic("K"),"=",6)))
			mtext(side=1,font=2,text=bquote(paste("true ",italic("K")," = ",3)),cex=1.2,padj=2.7)
		pie.plot(Q.file = "~/Dropbox/conStruct/sims/structure/datasets/simK3/simK3_estK.7.meanQ",
				 coords = sim.dataset$coords,
				 pop.vec = unlist(lapply(1:sim.dataset$N,function(n){rep(n,10)})),
				 layer.colors = laycols[[6]])
			mtext(side=3,font=2,text=bquote(paste(italic("K"),"=",7)))
dev.off()

