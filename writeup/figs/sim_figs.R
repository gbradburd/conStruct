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
pdf(file="~/Dropbox/conStruct/writeup/figs/sims/simK1_std.xval.pdf",width=6,height=6,pointsize=14)
	plot.xval.CIs(xval.CIs,K,simK=" (K = 1)")
	legend(x="bottomright",pch=19,col=c("blue","green"),legend=c("spatial","nonspatial"))
dev.off()

pdf(file="~/Dropbox/conStruct/writeup/figs/sims/simK1_std.xval_sp.pdf",width=6,height=6,pointsize=14)
	plot.xval.CIs(xval.CIs,K,simK=" (K = 1)",ylim=c(-10,0))
	legend(x="topright",pch=c(19,NA),lty=c(NA,1),lwd=c(NA,2),col=c(1,adjustcolor(1,0.8)),legend=c("mean","95% CI"))
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

#K2
setwd("~/Dropbox/conStruct/sims/cross_validation/K_2/x_validation")
n.reps <- 10
K <- 7
for(n in 1:n.reps){
	load(sprintf("simK2_rep%s_test.lnl.Robj",n))
	assign(paste0("tl",n),test.lnl)
}
x.vals <- lapply(1:n.reps,function(n){get(sprintf("tl%s",n))})
x.vals.std <- lapply(x.vals,standardize.xvals)
xval.CIs <- get.xval.CIs(x.vals.std,K)
pdf(file="~/Dropbox/conStruct/writeup/figs/sims/simK2_std.xval.pdf",width=6,height=6,pointsize=14)
	plot.xval.CIs(xval.CIs,K,simK=" (K = 2)")
	legend(x="bottomright",pch=19,col=c("blue","green"),legend=c("spatial","nonspatial"))
dev.off()

pdf(file="~/Dropbox/conStruct/writeup/figs/sims/simK2_std.xval_zoom.pdf",width=6,height=6,pointsize=14)
	plot.xval.CIs(xval.CIs,K,simK=" (K = 2)",k.range=c(2:7))
	legend(x="bottomright",pch=19,col=c("blue","green"),legend=c("spatial","nonspatial"))
dev.off()

pdf(file="~/Dropbox/conStruct/writeup/figs/sims/simK2_std.xval_sp.pdf",width=6,height=6,pointsize=14)
	plot.xval.CIs(xval.CIs,K,simK=" (K = 2)",k.range=c(2:7),ylim=c(-9,0))
	legend(x="bottomleft",pch=c(19,NA),lty=c(NA,1),lwd=c(NA,2),col=c(1,adjustcolor(1,0.8)),legend=c("mean","95% CI"))
dev.off()