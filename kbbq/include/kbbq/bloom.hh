#ifndef KBBQ_BLOOM_HH
#define KBBQ_BLOOM_HH
#include <cstdint>
#include <utility>
#include "yak/yak.h"
#include "yak/yak-priv.h" //provides yak_hash64()
#include "include/htslib/htslib/hts.h"
#include "minionrng/minion.h"

#define PREFIXBITS 10
//the number of prefix hashes is 1<<PREFIXBITS - 1

/* 
	A C++ wrapper around Heng Li's bloom filter in yak:
	https://github.com/lh3/yak
 */
//YAK_BLK_SHIFT= 9
//YAK_BLK_MASK = ((1<<(YAK_BLK_SHIFT)) - 1)= 11111111
//                                       3 =       11
//                                       4 =      100
//                                       7 =      111
//                                       9 =     1001
//                   (1<<(YAK_BLK_SHIFT-3) =   100000
//                                      31 =    11111
// BF_SHIFT is the size of the bloom filter (2 ^ BF_SHIFT bits.)

namespace bloom{

typedef std::array<Bloom,((1<<PREFIXBITS)-1)> bloomary_t;

static inline uint64_t hash(uint64_t key, uint64_t mask){yak_hash64(key, mask);}

//hash each kmer in seq and return a vector of hashes.
//the returned vector is not necessarily of length seq.length() - k + 1
//since non-ATGC characters will not be hashable.
static std::vector<uint64_t> hash_seq(std::string seq, int k);

//subsample the hashes and insert to bloom filter b.
// static void subsample_and_insert(Bloom& b, std::vector<uint64_t> hashes, double alpha, minion::Random& rng);

//subsample the hashes and insert into the proper bloom filter based on the prefix.
//htsiter::KmerSubsampler is the prefferred way to do this.
static void subsample_and_insert(bloomary_t bfs, std::vector<uint64_t> hashes, double alpha);

std::array<std::vector<int>,2> overlapping_kmers_in_bf(std::string seq, bloomary_t b);

//return the total number of kmers in b
int nkmers_in_bf(std::string seq, bloomary_t b, int k);

//return the INCLUSIVE indices bounding the largest stretch of trusted sequence
//if the first value is -1, there are no trusted kmers.
//if the 2nd value is -1 (== std::string::npos), until the end of the string is trusted.
//thus the whole string being trusted looks like {0, std::string::npos}
//while no part of the string being trusted looks like {std::string::npos, std::string::npos}
std::array<size_t,2> find_longest_trusted_seq(std::string seq, bloomary_t b, int k);

//find the longest possible fix for the kmer at position (k-1) until the end
//return the best character (multiple in case of a tie) and the length of the fix.
//if the length of the fix is 0, no fix was found and correction should end.
std::pair<std::vector<char>, int> find_longest_fix(std::string seq, bloomary_t& t, int k)

class Bloom
{
public:
	// the constructor approximates yak_bf_init()
	Bloom(): Bloom(22, 4){} //approx 512MB
	// across the 2^10 filters (2^PREFIXBITS), use 2^22 bits each. 2^22 * 2^10 = 2^32 bits = 512 MB
	Bloom(int nshift): Bloom(nshift, 4){}
	Bloom(int nshift, int nhashes);
	// this is similar to yak_bf_destroy()
	~Bloom();
	int nshift; //size of hash; 9 <= nshift <= 55
	int nhashes; //number of hash functions; 4 seems OK
	uint8_t* bloom;
	//this is similar to yak_bf_insert()
	int insert(unsigned long long hash); //try to insert the hash. return number of hashes that hit
	int query_n(unsigned long long hash); //return the number of hashes that hit but don't insert.
	inline bool query(unsigned long long hash);
	double fprate(unsigned long long n); //fprate given approximate number of times bf was loaded.
	static int optimal_nhashes(int shift, int n); //calculate the optimal nhashes for a size and number of times loaded.
};

}
#endif