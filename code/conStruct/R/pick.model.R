pick.stan.model <- function(spatial,k){
	if(k == 1 & !spatial){
		return(stanmodels$oneK)
	} else if(k == 1 & spatial){
		return(stanmodels$space_oneK)
	} else if(k > 1 & !spatial){
		return(stanmodels$multiK)
	} else if(k > 1 & spatial){
		return(stanmodels$space_multiK)
	}
}
