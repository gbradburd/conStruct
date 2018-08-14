functions {
	matrix spCov(int N, real a0, real aD, real a2, matrix D, real phi){
		matrix[N,N] cov;
		for(i in 1:N){
			for(j in i:N){
				cov[i,j] = a0 * exp( -(aD* D[i,j])^a2) + phi;
				cov[j,i] = cov[i,j];
			}
		}
		return cov;
	}
	matrix admixed_covariance(int N, int K, vector alpha0, vector alphaD, vector alpha2, matrix geoDist, matrix w_mat, vector nugget, vector phi, real gamma) {
		matrix[N,N] parCov;
		matrix[N,N] Nug_mat;
		parCov = rep_matrix(0,N,N);
		Nug_mat = diag_matrix(nugget);
		for(k in 1:K){
			parCov = parCov + tcrossprod(to_matrix(w_mat[,k])) .* spCov(N,alpha0[k],alphaD[k],alpha2[k],geoDist,phi[k]);
		}
		parCov = gamma + parCov + Nug_mat;
		return parCov;	
	}
	matrix make_w_matrix(int N, int K, vector[] w){
		matrix[N,K] w_mat;
		for(i in 1:N){
			w_mat[i] = to_row_vector(w[i]);
		}
		return w_mat;
	}
}
data {
	int<lower=1> K;		  				// number of layers
	int<lower=2> N; 	  				// number of samples
	int<lower=N+1> L;	    			// number of loci
	matrix[N,N] obsCov; 				// observed projected covariance
	matrix[N, N] geoDist; 				// matrix of pairwise geographic distance
	real varMeanFreqs;
	real<lower=0,upper=1> temp;			// temperature parameter for estimating marginal likelihood	
}
transformed data {
	matrix[N,N] LobsCov;				// n.loci multiplied by the sample covariance
	vector[K] dirConPar;
	LobsCov  = L * obsCov;
	dirConPar = rep_vector(0.1,K);
}
parameters {
	vector<lower=0>[K] alpha0;								// sill of the parametric covariance in layer k
	vector<lower=0>[K] alphaD;								// effect of geographic distance in the parametric covariance in layer k
	vector<lower=0, upper=2>[K]  alpha2;					// exponential slope parameter in the parametric covariance in layer k
	positive_ordered[K] phi;									// shared drift effect in layer k
  	vector<lower=0>[N] nugget; 								// sample-specific variance (allele sampling error + sample-specific drift)
	simplex[K]    w[N];    									// every sample (N in total) has a K simplex (i.e. K layers)
	real<lower=0> gamma;
}
transformed parameters {
	matrix[N,N] parCov;					// this specifies the parametric, admixed covariance matrix
	matrix[N,K] w_mat;
	w_mat = make_w_matrix(N,K,w);
	parCov = admixed_covariance(N, K, alpha0, alphaD, alpha2, geoDist, w_mat, nugget, phi, gamma);
}
model {
	alpha0 ~ normal(0,1);										// prior on alpha0
	alphaD ~ normal(0,1);										// prior on alphaD
	alpha2 ~ uniform(0,2);										// prior on alpha2
	nugget ~ normal(0,1);										// prior on nugget
	phi ~ normal(0,1);
	gamma ~ normal(varMeanFreqs,0.5);
	for(i in 1:N) w[i] ~ dirichlet(dirConPar);		// prior on admixture proportions
	LobsCov ~ wishart(L,parCov);					// likelihood function
}
