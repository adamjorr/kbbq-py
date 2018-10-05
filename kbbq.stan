data {
	int<lower = 1> L; //number of reads
	int<lower = 1> seqlen; //length of sequence
	int<lower = 1, upper = seqlen> ksize; //kmer size
	real<lower = 0> alpha;
	real<lower = 0> beta;
	int a[L, seqlen - ksize + 1]; //observed kmer abundances
}
transformed data{
	int nkmers; //number of kmers in each read
	nkmers = seqlen - ksize + 1;
}
parameters {
	real<lower = 0, upper = 1> e[L, seqlen]; //base error probability
	//p(kerr) could be NB or BINOM ?? P(nerr >= 1 given p(err), multiple trials) BETA BINOM???
	positive_ordered[2] lambda;
}
transformed parameters {
	real<upper = 0> kerr[L, nkmers]; //log probability kmer has 0 errors
	for (l in 1:L){
		for(k in 1:nkmers)
			kerr[l,k] = sum(log1m(segment(e[l],k,ksize))); //start at k, take slice of size ksize
	}
}
model {
	for (l in 1:L){
		for (n in 1:nkmers){
			target += log_sum_exp(kerr[l,n] + poisson_lpmf(a[l,n] | lambda[1]),
									log1m_exp(kerr[l,n]) + poisson_lpmf(a[l,n] | lambda[2]) );
		}
	} // abundances are drawn from this mixture
}
