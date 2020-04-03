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

kmer_cache_t find_trusted_kmers(HTSFile* file, bloom::bloomary_t& sampled, std::vector<int> thresholds, uint64_t chunksize = -1){
	uint64_t counted;
	kmer_cache_t ret;
	uint64_t x[2];
	uint64_t mask = (1ULL<<k*2) - 1;
	uint64_t shift = (k-1) * 2;
	int i, l;
	//the order here matters since we don't want to advance the iterator if we're chunked out
	while(counted++ < chunksize && file->next() >= 0){
		readutils::CReadData read = file->get();
		read.infer_errors(sampled, thresholds)
		int n_trusted = 0;
		for(i = l = 0, x[0] = x[1] = 0; i < read.seq.length(); ++i){
			if(!read->errors[i]){
				n_trusted++;
			}
			if(i >= k && !read.errors[i-k]){
				n_trusted--;
			}
			int c = seq_nt4_table[read->seq[i]];
			if (c < 4){
				x[0] = (x[0] << 2 | c) & mask;                  // forward strand
				x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;  // reverse strand
				if (++l >= k) { // we find a k-mer
					uint64_t y = x[0] < x[1]? x[0] : x[1];
					if(n_trusted == k){
						ret[y&((1<<PREFIXBITS)-1)].push_back(y);
					}
				}
			} else l = 0, x[0] = x[1] = 0; // if there is an "N", restart
		}
	}
	return ret;
}

covariateutils::CCovariateData get_covariatedata(HTSFile* file, bloom::bloomary_t& trusted, int k){
	CCovariateData data;
	while(file->next() >= 0){
		readutils::CReadData read = file->get();
		read.get_errors(trusted, k);
		data.consume_read(read);
	}
	return data;
}

}