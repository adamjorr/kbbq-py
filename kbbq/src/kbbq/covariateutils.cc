#include "covariateutils.hh"
#include <array>
#include "include/htslib/htslib/sam.h"

namespace covariateutils{

	void CCovariate::increment(size_t idx, const covariate_t& value){
		this[idx][0] += value[0]
		this[idx][1] += value[1]
	}

	void CRGCovariate::consume_read(const CReadData& read){
		int rg = read.get_rg_int();
		this->resize(rg);
		std::vector<bool> nse = read.not_skipped_errors();
		for(size_t i = 0; i < read.skips.length(); ++i){
			this->increment(std::array<int,2>((int)nse[i], (int)read.skips[i]));
		}
	}

	void CQCovariate::consume_read(const CReadData& read){
		int rg = read.get_rg_int();
		int q;
		this->resize(rg);
		for(size_t i = 0; i < read.skips.length(); ++i){
			if(!read.skips[i]){
				q = read.qual[i]
				this[rg].resize(q)
				this[rg][q].increment(std::array<int,2>((int)read.errors[i], (int)read.skips[i]));
			}
		}
	}

	void CCycleCovariate::consume_read(const CReadData& read){
		int rg = read.get_rg_int();
		int q;
		int cycle;
		this->resize(rg);
		size_t readlen = read.skips.length();
		for(size_t i = 0; i < readlen; ++i){
			if(!read.skips[i]){
				q = read.qual[i]
				this[rg].resize(q)
				this[rg][q][read.second].resize(i);
				this[rg][q][read.second][i].increment(std::array<int,2>((int)read.errors[i], 1));
			}
		}
	}

	//make sure first and second are 2 bit encoded nucleotides.
	//ensure there's no N!!!!
	inline int8_t dinuc_to_int(char first, char second){
		return 15 & ((first << 2) | second); // 1111 & (xx00 | 00xx)
	}

	inline std::array<char, 2> int_to_dinuc(int dinuc){
		std::array<char,2> x = {dinuc >> 2, dinuc | 3} ;
		return x;
	}

	void CDinucCovariate::consume_read(const CReadData& read, int minscore = 6){
		int rg = read.get_rg_int();
		this->resize(rg);
		int q;
		int8_t dinuc;
		for(size_t i = 1; i < read.seq.length(); i++){
			q = read.qual[i];
			if(!read.skips[i] && read.seq[i-1] != 'N' && read.seq[i] != 'N' && q >= minscore){
				this[rg].resize(q);
				char first = seq_nt16_int[seq_nt16_table[read.seq[i-1]]];
				char second = seq_nt16_int[seq_nt16_table[read.seq[i]]];
				this[rg][q][dinuc_to_int(first, second)].increment(std::array<int,2>((int)read.errors[i], 1))
			}
		}
		// seq_nt16_table[256]: char -> 4 bit encoded (1/2/4/8)
		// seq_nt16_str[]: 4 bit -> char
		// seq_nt16_int[]: 4 bit -> 2 bits (0/1/2/3)
	}

	void CCovariateData::consume_read(const CReadData& read, int minscore){
		rgcov.consume_read(read);
		qcov.consume_read(read);
		cycov.consume_read(read);
		dicov.consume_read(read, minscore);
	}





}