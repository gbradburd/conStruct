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