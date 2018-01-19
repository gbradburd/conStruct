for(k in 1:7){
	infile <- "poplar.bed"
	call <- paste0("admixture --cv=50 ",infile ," ",k)
	system(call)
}