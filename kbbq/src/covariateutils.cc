#include "covariateutils.hh"
#include <array>

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
				this[rg][q][read.second][i].increment(std::array<int,2>((int)read.errors[i], (int)read.skips[i]));
			}
		}
	}

	void CDinucCovariate::consume_read(const CReadData& read){
		//TODO
	}

	void CCovariateData::consume_read(const CReadData& read){
		//TODO
	}





}