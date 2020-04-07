#ifndef RECALIBRATEUTILS_H
#define RECALIBRATEUTILS_H

#include <vector>
#include <string>
#include <cmath>
#include "readutils.hh"
#include "covariateutils.hh"
#include "bloom.hh"
#include "htsiter.hh"

namespace recalibrateutils{

typedef meanq_t std::vector<int>;
typedef globaldq_t std::vector<int>;
typedef qscoredq_t std::vector<std::vector<int>>;
typedef cycledq_t std::vector<std::vector<std::vector<int>>>;
typedef dinucdq_t std::vector<std::vector<std::vector<int>>>;

struct dq_t
{
	meanq_t meanq;
	globaldq_t globaldq;
	dqscoredq_t qscoredq;
	cycledq_t cycledq;
	dinucdq_t dinucdq;
};

typedef std::array<std::vector<uint64_t>, 1<<PREFIXBITS-1> kmer_cache_t;

//subsample kmers, hash them, and add them to a prefix tree
//read chunksize kmers at once. if chunksize = -1, read up to MAXINT
kmer_cache_t subsample_kmers(KmerSubsampler& s, uint64_t chunksize = -1);

//read all kmers from the subsampler and load them into a cache.
kmer_cache_t read_all_kmers(KmerSubsampler& s, uint64_t chunksize = -1);

//add kmers in a given cache to the given assembly of bloom filters.
void add_kmers_to_bloom(kmer_cache_t& kmers, bloom::bloomary_t& filters);

//get some reads from a file, whether a kmer is trusted and put it in a cache.
//this can probably be parallelized if needed because the reads are independent
kmer_cache_t find_trusted_kmers(HTSFile* file, bloom::bloomary_t& sampled, std::vector<int> thresholds, uint64_t chunksize = -1);

inline long double q_to_p(int q){return std::powl(10.0l, -((long double)q / 10.0l))}
inline int p_to_q(long double p, int maxscore = 42){return p > 0 ? (int)(-10 * std::log10(p)) : maxscore}

//get covariate data using the trusted kmers
//this can be parallelized easily since each read is independent
covariateutils::CCovariateData get_covariatedata(HTSFile* file, bloom::bloomary_t& trusted);

//recalibrate all reads given the CovariateData
void recalibrate_and_write(HTSFile* in, dq_t dqs, std::string outfn);

}

#endif