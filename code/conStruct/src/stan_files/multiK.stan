functions {
	matrix admixed_covariance(int N, int K, matrix w_mat, vector nugget, vector phi, real gamma) {
		matrix[N,N] parCov;
		matrix[N,N] Nug_mat;
		parCov = rep_matrix(0,N,N);
		Nug_mat = diag_matrix(nugget);
		for(k in 1:K){
			parCov = parCov + tcrossprod(to_matrix(w_mat[,k])) * phi[k];
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
	positive_ordered[K] phi;				// shared drift effect in layer k
	real<lower=0> gamma;				// covariance between all layers
  	vector<lower=0>[N] nugget; 			// sample-specific variance (allele sampling error + sample-specific drift)
	simplex[K]    w[N];    				// every sample (N in total) has a K simplex (i.e. K layers)
}
transformed parameters {
	matrix[N,N] parCov;					// this specifies the parametric, admixed covariance matrix
	matrix[N,K] w_mat;
	w_mat = make_w_matrix(N,K,w);
	parCov = admixed_covariance(N, K, w_mat, nugget, phi, gamma);
}
model {
	nugget ~ normal(0,1);										// prior on nugget
	phi ~ normal(0,1);
	gamma ~ normal(varMeanFreqs,0.5);
	for(i in 1:N) w[i] ~ dirichlet(dirConPar);				    // prior on admixture proportions
	LobsCov ~ wishart(L,parCov);						// likelihood function
}

