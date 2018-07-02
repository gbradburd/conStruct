for(k in 1:7){
	infile <- "simK1.bed"
	call <- paste0("admixture --cv=50 ",infile ," ",k)
	system(call)
}