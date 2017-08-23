collapse.ind.Q <- function(Q,pop.vec){
	n.pops <- length(unique(pop.vec))
	admix.props <- matrix(NA,n.pops,ncol(Q))
	for(i in 1:n.pops){
		admix.props[i,] <- colMeans(Q[which(pop.vec==i),])
	}
	return(admix.props)
}

pie.plot <- function(Q.file,coords,pop.vec,cluster.colors){
	Q <- read.table(Q.file,stringsAsFactors=FALSE)
	admix.props <- collapse.ind.Q(Q,pop.vec)
	K <- ncol(admix.props)
	N <- length(unique(pop.vec))
	cluster.names <- paste0("cluster_",1:K)
	sample.names <- paste0("sample_", 1:N)
	color.tab <- caroline::nv(c(cluster.colors[1:K]), cluster.names)
	pie.list <- lapply(1:N, function(i) {
	    caroline::nv(admix.props[i, ], cluster.names)
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

load("~/Dropbox/conStruct/sims/cross_validation/K_1/sim.dataset.Robj")

pdf(file="~/Dropbox/conStruct/writeup/figs/fastStr/fastStr_simK1_pies.pdf",width=12,height=4.5,pointsize=13)
	par(mfrow=c(1,3))
		pie.plot(Q.file = "~/Dropbox/conStruct/sims/structure/datasets/simK1/simK1_estK2.2.meanQ",
				 coords = sim.dataset$coords,
				 pop.vec = unlist(lapply(1:sim.dataset$N,function(n){rep(n,10)})),
				 cluster.colors = c(4,2))
			mtext(side=3,font=2,text="K=2")
		pie.plot(Q.file = "~/Dropbox/conStruct/sims/structure/datasets/simK1/simK1_estK3.3.meanQ",
				 coords = sim.dataset$coords,
				 pop.vec = unlist(lapply(1:sim.dataset$N,function(n){rep(n,10)})),
				 cluster.colors = c(2,"goldenrod1",4))
			mtext(side=3,font=2,text="K=3")
			mtext(side=1,font=2,text="true K = 1",cex=1.2,padj=2.7)
		pie.plot(Q.file = "~/Dropbox/conStruct/sims/structure/datasets/simK1/simK1_estK4.4.meanQ",
				 coords = sim.dataset$coords,
				 pop.vec = unlist(lapply(1:sim.dataset$N,function(n){rep(n,10)})),
				 cluster.colors = c("goldenrod1",4,2,"forestgreen"))
			mtext(side=3,font=2,text="K=4")
dev.off()


load("~/Dropbox/conStruct/sims/cross_validation/K_2/sim.dataset.Robj")

pdf(file="~/Dropbox/conStruct/writeup/figs/fastStr/fastStr_simK2_pies.pdf",width=12,height=4.5,pointsize=13)
	par(mfrow=c(1,3))
		pie.plot(Q.file = "~/Dropbox/conStruct/sims/structure/datasets/simK2/simK2_est2.2.meanQ",
				 coords = sim.dataset$coords,
				 pop.vec = unlist(lapply(1:sim.dataset$N,function(n){rep(n,10)})),
				 cluster.colors = c(4,2))
			mtext(side=3,font=2,text="K=2")
		pie.plot(Q.file = "~/Dropbox/conStruct/sims/structure/datasets/simK2/simK2_estK3.3.meanQ",
				 coords = sim.dataset$coords,
				 pop.vec = unlist(lapply(1:sim.dataset$N,function(n){rep(n,10)})),
				 cluster.colors = c(4,2,"goldenrod1"))
			mtext(side=3,font=2,text="K=3")
			mtext(side=1,font=2,text="true K = 2",cex=1.2,padj=2.7)
		pie.plot(Q.file = "~/Dropbox/conStruct/sims/structure/datasets/simK2/simK2_estK4.4.meanQ",
				 coords = sim.dataset$coords,
				 pop.vec = unlist(lapply(1:sim.dataset$N,function(n){rep(n,10)})),
				 cluster.colors = c(2,4,"goldenrod1","forestgreen"))
			mtext(side=3,font=2,text="K=4")
dev.off()
			
load("~/Dropbox/conStruct/sims/cross_validation/K_3/sim.dataset.Robj")
pdf(file="~/Dropbox/conStruct/writeup/figs/fastStr/fastStr_simK3_pies.pdf",width=12,height=4.5,pointsize=13)
	par(mfrow=c(1,3))
		pie.plot(Q.file = "~/Dropbox/conStruct/sims/structure/datasets/simK3/simK3_estK2.2.meanQ",
				 coords = sim.dataset$coords,
				 pop.vec = unlist(lapply(1:sim.dataset$N,function(n){rep(n,10)})),
				 cluster.colors = c(4,2))
			mtext(side=3,font=2,text="K=2")
		pie.plot(Q.file = "~/Dropbox/conStruct/sims/structure/datasets/simK3/simK3_estK3.3.meanQ",
				 coords = sim.dataset$coords,
				 pop.vec = unlist(lapply(1:sim.dataset$N,function(n){rep(n,10)})),
				 cluster.colors = c(4,"goldenrod1",2))
			mtext(side=3,font=2,text="K=3")
			mtext(side=1,font=2,text="true K = 3",cex=1.2,padj=2.7)
		pie.plot(Q.file = "~/Dropbox/conStruct/sims/structure/datasets/simK3/simK3_estK4.4.meanQ",
				 coords = sim.dataset$coords,
				 pop.vec = unlist(lapply(1:sim.dataset$N,function(n){rep(n,10)})),
				 cluster.colors = c(2,4,"goldenrod1","forestgreen"))
			mtext(side=3,font=2,text="K=4")
dev.off()
