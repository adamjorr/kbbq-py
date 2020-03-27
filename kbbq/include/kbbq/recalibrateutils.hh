#ifndef RECALIBRATEUTILS_H
#define RECALIBRATEUTILS_H

#include <vector>
#include <string>
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

std::vector<int> recalibrate_read(CReadData read, dq_t dqs, int minscore = 6);

typedef std::array<std::vector<uint64_t>, 1<<PREFIXBITS-1> kmer_cache_t;

//subsample kmers, hash them, and add them to a suffix tree
//read chunksize kmers at once. if chunksize = 0, read them all
kmer_cache_t subsample_kmers(KmerSubsampler& s, int chunksize = 0);

//add kmers in a given cache to the given assembly of bloom filters.
void add_kmers_to_bloom(kmer_cache_t kmers, bloom::bloomary_t& filters);

}

#endif