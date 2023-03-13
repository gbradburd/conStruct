functions {
	matrix spCov(int N, real a0, real aD, real a2, matrix D, vector nugget, real gamma) {
		matrix[N,N] parCov;
		matrix[N,N] Nug_mat;
		parCov = rep_matrix(0,N,N);
		Nug_mat = diag_matrix(nugget);
		for(i in 1:N){
			for(j in i:N){
				parCov[i,j] = a0 * exp( -(aD * D[i,j])^a2);
				parCov[j,i] = parCov[i,j];
			}
		}
		parCov += gamma + Nug_mat;
		return parCov;	
	}
}
data {
	int<lower=1> K;		  				// number of layers
	int<lower=2> N; 	  				// number of samples
	int<lower=N+1> L;	    			// number of loci
	matrix[N,N] obsCov; 				// observed projected covariance
	matrix[N, N] geoDist; 				// matrix of pairwise geographic distance 
	real varMeanFreqs;
}
transformed data {
	matrix[N,N] LobsCov;				// n.loci multiplied by the sample covariance
	LobsCov  = L * obsCov;
}
parameters {
	real<lower=0> alpha0;								// sill of the parametric covariance in layer k
	real<lower=0> alphaD;								// effect of geographic distance in the parametric covariance in layer k
	real<lower=0, upper=2>  alpha2;					// exponential slope parameter in the parametric covariance in layer k
	real<lower=0> gamma;				// covariance between all layers
  	vector<lower=0>[N] nugget; 								// sample-specific variance (allele sampling error + sample-specific drift)
}
transformed parameters {
	matrix[N,N] parCov;					// this specifies the parametric, admixed covariance matrix
	parCov = spCov(N, alpha0, alphaD, alpha2, geoDist, nugget, gamma);
}
model {
	alpha0 ~ normal(0,1);										// prior on alpha0
	alphaD ~ normal(0,1);										// prior on alphaD
	alpha2 ~ uniform(0,2);										// prior on alpha2
	nugget ~ normal(0,1);										// prior on nugget
	gamma ~ normal(varMeanFreqs,0.5);							// prior on global covariance
	LobsCov ~ wishart(L,parCov);								// likelihood function
}
