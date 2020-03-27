#include "recalibrateutils.hh"
#include "htsiter.hh"

namespace recalibrateutils{

//if the number of reads doesn't fit in a uint64_t call this twice :)
kmer_cache_t subsample_kmers(KmerSubsampler& s, uint64_t chunksize = -1){
	uint64_t kmer = 0;
	uint64_t counted;
	kmer_cache_t ret;
	//the order here matters since we don't want to advance the iterator if we're chunked out
	while(counted++ < chunksize && ((kmer = s.next_kmer()) != 0 || !s.readseq.empty())){
		prefix = kmer & ((1<<PREFIXBITS)-1);
		ret[prefix].push_back(kmer);
	}
	return ret;
}

//this is a good target for multithreading ;)
void add_kmers_to_bloom(kmer_cache_t& kmers, bloom::bloomary_t& filters){
	for(int i = 0; i < (1<<PREFIXBITS); ++i){
		for(uint64_t kmer : kmers[i]){
			filters[i].insert(kmer >> PREFIXBITS); //remove the suffix and put into the filter
		}
	}
}


}