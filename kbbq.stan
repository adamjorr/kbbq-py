data {
	int<lower = 1> L; //number of reads
	int<lower = 1> seqlen; //length of sequence
	int<lower = 1, upper = seqlen> ksize; //kmer size
	int a[L, seqlen - ksize + 1]; //observed kmer abundances
	int<lower = 1, upper = 43> q[L, seqlen]; //quality scores
}
transformed data{
	int nkmers; //number of kmers in each read
	real<lower = 0, upper = 1> initial_e[L, seqlen]; //initial probability of error from q score
	int<lower = 1> maxq; //maximum quality score in this data

	nkmers = seqlen - ksize + 1;
	for (l in 1:L){
		for (s in 1:seqlen){
			initial_e[l,s] = pow(10, -q[l,s] / 10.0);
		}
	}
	maxq = max(to_array_1d(q));
}
parameters {
	real<lower = 0, upper = 1> e[L, seqlen]; //base error probability
	//p(kerr) could be NB or BINOM ?? P(nerr >= 1 given p(err), multiple trials) BETA BINOM???
	positive_ordered[2] lambda; //parameters for mixed abundance distribution

	//we reparameterize alpha and beta into phi = alpha / (alpha + beta) and gamma = alpha + beta
	vector<lower = 0, upper = 1>[maxq] phi; //alpha[q] is the alpha for bases with quality score q
	vector<lower = 0.1>[maxq] gamma; //beta[q] is the beta for bases with quality score q
}
transformed parameters {
	real<upper = 0> kerr[L, nkmers]; //log probability kmer has 1 or more errors
	vector<lower = 0>[maxq] alpha;
	vector<lower = 0>[maxq] beta;

	for (l in 1:L){
		for(k in 1:nkmers)
			kerr[l,k] = log1m_exp(sum(log1m(segment(e[l],k,ksize)))); //start at k, take slice of size ksize
	}
	alpha = phi .* gamma;
	beta = gamma ./ (1 - phi);
}
model {
	phi ~ beta(2,2);
	gamma ~ gamma(2,1.0/100.0);
	lambda[1] ~ gamma(2,1);
	lambda[2] ~ gamma(2,1.0/500.0);
	for (l in 1:L){
		for (i in 1:seqlen){
			e[l,i] ~ beta(alpha[q[l,i]], beta[q[l,i]]);
		}
		for (n in 1:nkmers){
			target += log_sum_exp(kerr[l,n] + poisson_lpmf(a[l,n] | lambda[1]),
									log1m_exp(kerr[l,n]) + poisson_lpmf(a[l,n] | lambda[2]) );
		}
	} // abundances are drawn from this mixture
}

