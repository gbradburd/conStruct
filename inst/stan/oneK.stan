functions {
	matrix Cov(int N, vector nugget, real gamma) {
		matrix[N,N] parCov;
		matrix[N,N] Nug_mat;
		parCov = rep_matrix(gamma,N,N);
		Nug_mat = diag_matrix(nugget);
		parCov += Nug_mat;
		return parCov;	
	}
}
data {
	int<lower=1> K;		  				// number of layers
	int<lower=2> N; 	  				// number of samples
	int<lower=N+1> L;	    			// number of loci
	matrix[N,N] obsCov; 				// observed projected covariance
	real varMeanFreqs;					// variance in mean frequencies
}
transformed data {
	matrix[N,N] LobsCov;				// n.loci multiplied by the sample covariance
	LobsCov  = L * obsCov;
}
parameters {
	real<lower=0> gamma;				// covariance between all layers
  	vector<lower=0>[N] nugget; 			// sample-specific variance (allele sampling error + sample-specific drift)
}
transformed parameters {
	matrix[N,N] parCov;					// this specifies the parametric, admixed covariance matrix
	parCov = Cov(N, nugget, gamma);
}
model {
	nugget ~ normal(0,1);				// prior on nugget
	gamma ~ normal(varMeanFreqs,0.5);	// prior on gamma
	LobsCov ~ wishart(L,parCov);		// likelihood function
}
