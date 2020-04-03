#include <cstdlib>
#include <cstring>
#include <cmath>
#include <random>
#include "bloom.hh"
#include "minionrng/minion.hpp"

namespace bloom
{
	static std::vector<uint64_t> hash_seq(std::string seq, int k){
		int i, l;
		uint64_t x[2];
		uint64_t mask = (1ULL<<k*2) - 1;
		uint64_t shift = (k-1) * 2;
		std::vector<uint64_t> ret;
		ret.reserve(seq.length() - k + 1);
		for (i = l = 0, x[0] = x[1] = 0; i < seq.length(); ++i) {
			int c = seq_nt4_table[seq[i]];
			if (c < 4) { // not an "N" base
				x[0] = (x[0] << 2 | c) & mask;                  // forward strand
				x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;  // reverse strand
				if (++l >= k) { // we find a k-mer
					uint64_t y = x[0] < x[1]? x[0] : x[1];
					ret.push_back(yak_hash64(y, mask));
				}
			} else l = 0, x[0] = x[1] = 0; // if there is an "N", restart
		}
	}

	static void subsample_and_insert(bloomary_t bfs, std::vector<uint64_t> hashes, double alpha, minion::Random& rng){
		int mask = (1<<PREFIXBITS)-1;
		std::bernoulli_distribution dist(alpha);
		for(uint64_t h: hashes){
			if(dist(rng)){
				bfs[h&mask].insert(h >> PREFIXBITS);
			}
		}
	}

	std::array<std::vector<int>,2> overlapping_kmers_in_bf(std::string seq, bloomary_t b, int k = 31){
		int i, l; //i is the character index, l is the length of the current stretch of non-N bases
		uint64_t x[2];
		uint64_t mask = (1ULL<<k*2) - 1;
		uint64_t shift = (k-1) * 2;
		std::vector<bool> kmers; //kmer indices are in l space
		std::vector<int> inbf(seq.length()); //inbf and possible indices are in i space.
		std::vector<int> possible(seq.length());
		int n_in;
		int n_out;
		kmers.reserve(seq.length() - k + 1); //whether the kmer ENDING in that position is in the BF
		for (i = l = 0, x[0] = x[1] = 0; i < seq.length(); ++i) {
			int c = seq_nt4_table[seq[i]];
			if (c < 4) { // not an "N" base
				x[0] = (x[0] << 2 | c) & mask;                  // forward strand
				x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;  // reverse strand
				if (++l >= k) { // we find a k-mer
					// n_possible = n_possible + 1 <= k ? n_possible + 1 : k;
					// possible[i-k] = n_possible;
					uint64_t y = x[0] < x[1]? x[0] : x[1];
					uint64_t h = yak_hash64(y, mask);
					bool k_in = b[h&((1<<PREFIXBITS)-1)].query(h>>PREFIXBITS);
					kmers.push_back(k_in);
					if(k_in){
						n_in++; //n_in += !!k_in; n_out += !!!k_in
					} else{
						n_out++;
					}
					//we have a long enough stretch of kmers that the old ones start to matter
					if(l >= 2 * k){ //l-k>=k;
						if(kmers[l-k-1]){
							n_in--;
						} else {
							n_out--;
						}
					}

					inbf[i -k + 1] = n_in;
					possible[i -k + 1] = n_in + n_out; //n_possible = n_in + n_out
				} else {
					kmers.push_back(false);
				}
			} else{
				l = 0, x[0] = x[1] = 0, n_possible = 0; // if there is an "N", restart
				kmers.clear();
			}
		}
		if(l >= k){ //we ended with a full kmer so we have to deal with the end positions
			for(; i < seq.length()+k; ++i){
				if(++l >= 2 * k){ //l-k>=k;
					if(kmers[l-k-1]){
						n_in--;
					} else {
						n_out--;
					}
				}
				inbf[i - k + 1] = n_in;
				possible[i - k + 1] = n_in + n_out;
			}
		}
		std::array<std::vector<int>,2> ret = {inbf, possible};
		return ret;
	}

	int nkmers_in_bf(std::string seq, bloomary_t b, int k){
		std::vector<uint64_t> hashes = hash_seq(seq, k);
		int c = 0;
		for(uint64_t h: hashes){
			if(b[h&((1<<PREFIXBITS)-1)].query(h>>PREFIXBITS)){
				c++;
			}
		}
		return c;
	}

	std::array<size_t, 2> find_longest_trusted_seq(std::string seq, bloomary_t& b, int k){
		std::vector<uint64_t> hashes = bloom::hash_seq(seq, int k);
		int i, l, n;
		int anchor_start, anchor_end, anchor_l, anchor_current;
		anchor_start = anchor_end = -1;
		for(i = l = n = 0; i < seq.length(); ++i){
			int c = seq_nt4_table[seq[i]];
			if(c < 4){
				if(++l > k){
					uint64_t h = hashes[n++];
					if(b[h&((1<<PREFIXBITS)-1)].query(h>>PREFIXBITS)){
						anchor_current++; //length of current stretch
					} else {
						if(anchor_current > anchor_l){
							anchor_l = anchor_current;
							anchor_end = i;
							anchor_start = i + 1 - k - anchor_l;
						}
						anchor_current = 0;
					}
				}
			} else {
				if(anchor_current > anchor_l){
					anchor_l = anchor_current;
					anchor_end = i-1;
					anchor_start = i - k - anchor_l;
				}
				anchor_current = 0;
				l = 0;
			}
		}
		if(anchor_current > anchor_l){
			anchor_l = anchor_current;
			anchor_end = std::string::npos;
			anchor_start = i - k - anchor_l;
		}
		return std::array<size_t,2>{{anchor_start, anchor_end}};
	}

	void increment_coded_kmer(uint64_t x[], int c, int k){
		uint64_t mask = (1ULL<<k*2) - 1;
		uint64_t shift = (k-1) * 2;
		if(c < 4){
			x[0] = (x[0] << 2 | c) & mask;                  // forward strand
			x[1] = x[1] >> 2 | (uint64_t)(3 - c) << shift;  // reverse strand
		}
	}

	std::pair<std::vector<char>, int> find_longest_fix(std::string seq, bloomary_t& t, int k){
		//this function is probably not very efficient.
		size_t i;
		int l;
		int best_l = 0;
		uint64_t x[2];
		std::vector<char> best_c;
		char original_c;
		size_t end = std::min(seq.length(), 2 * k - 1);
		original_c = seq[k-1];
		for(char d : {"A","C","G","T"} ){
			if( d != original_c ){
				seq[k-1] = d;
				// current_l = nkmers_in_bf(seq.substr(0, end), t, k); //this is not very efficient but that's ok
				for(i = l = 0, x[0] = x[1] = 0; i < end; ++i){
					int c = seq_nt4_table[seq[i]];
					if (c < 4) { // not an "N" base
						increment_coded_kmer(x, c, k);
						if (++l >= k) { // we find a k-mer
							uint64_t y = x[0] < x[1]? x[0] : x[1];
							uint64_t h = yak_hash64(y, mask);
							if(!t[h&((1<<PREFIXBITS)-1)].query(h>>PREFIXBITS)){
								break; //if it's not trusted, end here
							}
						}
					} else {
						break;
					}
				}
				if(l > best_l){
					best_l = l;
					best_c.clear() 
					best_c.push_back(d);
				} else if (l == best_l){
					best_c.push_back(d);
				}
			}
		}
		if(best_l == seq.length() && best_c.size() > 1){ //we ran out of kmers to try and we have a tie
			for(char d : {"A","C","G","T"} ){
				if( d != original_c ){
					seq = seq.substr(0, end); //original sequence length
					seq[k-1] = d;
					for(i = end - k, l = 0, x[0] = x[1] = 0; i < 2 * k - 1; ++i){
						int cur_seqlen = seq.length();
						for(char extra : {"A","C","G","T"}){		
							int c = i < seq.length() ? seq_nt4_table[seq[i]] : seq_nt4_table[extra];
							if (c < 4) { // not an "N" base
								increment_coded_kmer(x, c, k);
								if (++l >= k) { // we find a k-mer
									uint64_t y = x[0] < x[1]? x[0] : x[1];
									uint64_t h = yak_hash64(y, mask);
									if(i >= seq.length() && t[h&((1<<PREFIXBITS)-1)].query(h>>PREFIXBITS)){
										seq.push_back(extra)
										break; //if we have an extra and it's trusted, end here
									}
								}
							}
						}
						if(cur_seqlen == seq.length()){ //we didn't add an extra base
							break; //exit extension
						}
					}
					if(l > best_l){
						best_l = l;
						best_c.clear() 
						best_c.push_back(d);
					} else if (l == best_l){
						best_c.push_back(d);
					}
				}
			}
		}
		return std::make_pair(best_c, k - best_l);
	}

	Bloom::Bloom(int nshift, int nhashes): nshift(nshift), nhashes(nhashes) {
		//nshift + YAK_BLK_SHIFT should be less than 64 (nshift <= 55)
		//nshift should be greater than or equal to YAK_BLK_SHIFT (9)
		//thus 9 <= nshift <= 55
		void* aligned = 0;
		posix_memalign(&aligned, 1<<(YAK_BLK_SHIFT-3), 1ULL<<(nshift-3)); //from cstdlib
		this->bloom = aligned;
		bzero(this->bloom, 1ULL<<(nshift-3));
	}

	Bloom::~Bloom(){
		delete this->bloom;
	}

	int Bloom::insert(unsigned long long hash){
		int x = this->nshift - YAK_BLK_SHIFT; // the bloom filter size
		unsigned long long y = hash & ((1ULL<<x)-1); //last x bits of hash;
		int h1 = hash >> x & YAK_BLK_MASK; // 8 bits right before y;
		int h2 = hash >> this->nshift & YAK_BLK_MASK; //8 bits right before h1;
		uint8_t *p = &this->bloom[y << (YAK_BLK_SHIFT-3)];
		if((h2&31) == 0){
			h2 = (h2 + 1) & YAK_BLK_MASK;
		}
		int z = h1;
		int count = 0;
		for(int i = 0; i < this->nhashes; z = (z + h2) & YAK_BLK_MASK){
			uint8_t *q = &p[z>>3];
			uint8_t u = 1 << (z&7);
			count += !!(*q & u); // !! = double negate; makes 1 or 0
			q |= u;
			++i;
		}
		return count;
	}

	int Bloom::query_n(unsigned long long hash){
		int x = this->nshift - YAK_BLK_SHIFT; // the bloom filter size
		unsigned long long y = hash & ((1ULL<<x)-1); // fit the hash into the bloom filter;
		//discard the beginning (which determines which filter the hash goes into to begin with)
		int h1 = hash >> x & YAK_BLK_MASK; // this is the beginning part that's not in y;
		int h2 = hash >> this->nshift & YAK_BLK_MASK; //this is some middle part of the hash;
		uint8_t *p = &this->bloom[y << (YAK_BLK_SHIFT-3)];
		if((h2&31) == 0){
			h2 = (h2 + 1) & YAK_BLK_MASK;
		}
		int z = h1;
		int count = 0;
		for(int i = 0; i < this->nhashes; z = (z + h2) & YAK_BLK_MASK){
			uint8_t *q = &p[z>>3];
			uint8_t u = 1 << (z&7);
			count += !!(*q & u); // !! = double negate; makes 1 or 0
			++i;
		}
		return count;
	}

	inline bool Bloom::query(unsigned long long hash){
		return this->query_n(hash) == this->nhashes;
	}

	//fprate given approximate number of times bf was loaded
	double Bloom::fprate(unsigned long long n){
		int m = 1ULL<<(this->nshift); //number of bits in the filter
		int k = this -> nhashes;
		return pow(1 - exp(-1.0 * k * n / m), k);
	}

	//there are 2 ** shift bits in the filter.
	static int Bloom::optimal_nhashes(int shift, int n){
		return floor(1.0 * (2 ** shift) / n * log(2))
	}

	// bool Bloom::kmers_in(std::string seq, int k){
	// 	std::vector<uint64_t> hashes = hash_seq(seq, k);
	// 	std::vector<bool> in(hashes.size());
	// 	for(int i = 0; i < hashes.size(); i++){
	// 		in[i] = this->query(hashes[i]);
	// 	}
	// }

	/* //let's consider the case when nshift is minimal; or nshift = YAK_BLK_SHIFT (9)
	Bloom::insert(unsigned long long hash){
		int h1 = hash & YAK_BLK_MASK; //last 8 bits
		int h2 = hash >> 9 & YAK_BLK_MASK; //skip one then get the second to last 8 bits
		// this looks like:
		//... [h2] [h2] [h2] [h2] [h2] [h2] [h2] [h2] [] [h1] [h1] [h1] [h1] [h1] [h1] [h1] [h1]
		//                                                 z    z    z    z    z <- initial address in bf
		//                                                            chosen bit -> z    z    z     
		uint8_t *p = &this->bloom[0];
		if((h2&31) == 0){ //if last 5 bits of h2 == 0;
			h2 = (h2 + 1) & YAK_BLK_MASK; //now last 5 bits == 1
		}
		int z = h1; // z = last 8 bits of hash
		int count = 0; //initialize count
		for(int i = 0; i < this->nhashes; z = (z + h2) & YAK_BLK_MASK){ //for every hash fn;
			uint8_t *q = &p[z >> 3]; //q = pointer to some part of the bloom filter
			uint8_t u = 1 << (z&7); // 1 << (last 3 bits of z); the bit chosen to look at
			count += !!(*q & u); // increment count if it's there
			q |= u; //update filter
			++i; // go to next hash
			//now z = z + h2
		}


	} */
}






