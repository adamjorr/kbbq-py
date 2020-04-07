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

	//TODO: fix priors
	rgdq_t CRGCovariate::delta_q(meanq_t prior){
		rgdq_t dq(this->size());
		for(int i = 0; i < this->size(); ++i){ //i is rgs here
			int map_q = 0; //maximum a posteriori q
			int best_posterior = 0;
			for(int possible = 0; possible < MAXQ+1; possible++){
				int diff = std::abs(prior[i] - possible);
				long double prior_prob = normal_prior[diff];
				long double p = recalibrateutils::q_to_p(possible);
				long double loglike = log_binom_pmf(*this[i][0], *this[i][1], p);
				long double posterior = prior_prob + loglike;
				if(posterior > best_posterior){
					map_q = possible;
					best_posterior = posterior;
				}
			}
			dq[i] = map_q - prior[i];
		}
		return dq;
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

	qscoredq_t CQCovariate::delta_q(prior1_t prior){
		qscoredq_t dq(this->size());
		for(int i = 0; i < this->size(); ++i){ //i is rgs here
			dq[i] = std::vector(*this[i].size());
			for(int j = 0; j < *this[i].size(); ++j){ //j is q's here
				int map_q = 0; //maximum a posteriori q
				int best_posterior = 0;
				for(int possible = 0; possible < MAXQ+1; possible++){
					int diff = std::abs(prior[i] - possible);
					long double prior_prob = normal_prior[diff];
					long double p = recalibrateutils::q_to_p(possible);
					long double loglike = log_binom_pmf(*this[i][j][0], *this[i][j][1], p);
					long double posterior = prior_prob + loglike;
					if(posterior > best_posterior){
						map_q = possible;
						best_posterior = posterior;
					}
				}
				dq[i][j] = map_q - prior[i];
			}
		}
		return dq;
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

	cycledq_t CCycleCovariate::delta_q(prior2_t prior){
		cycledq_t dq(this->size()); //rg -> q -> fwd(0)/rev(1) -> cycle -> values
		for(int i = 0; i < this->size(); ++i){ //i is rgs here
			dq[i] = std::vector(*this[i].size());
			for(int j = 0; j < *this[i].size(); ++j){ //j is q's here
				for(int k = 0; k < 2; ++k){ //fwd/rev
					dq[i][j][k] = std::vector(*this[i][j][k].size());
					for(l = 0; l < *this[i][j][k].size()){ //cycle value
						int map_q = 0; //maximum a posteriori q
						int best_posterior = 0;
						for(int possible = 0; possible < MAXQ+1; possible++){
							int diff = std::abs(prior[i][j] - possible);
							long double prior_prob = normal_prior[diff];
							long double p = recalibrateutils::q_to_p(possible);
							long double loglike = log_binom_pmf(*this[i][j][k][l][0], *this[i][j][k][l][1], p);
							long double posterior = prior_prob + loglike;
							if(posterior > best_posterior){
								map_q = possible;
								best_posterior = posterior;
							}
						}
						dq[i][j][k][l] = map_q - prior[i][j];
					}
				}
			}
		}
		return dq;
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

	dinucdq_t CDinucCovariate::delta_q(prior2_t prior){
		dinucdq_t dq(this->size()); //rg -> q -> fwd(0)/rev(1) -> cycle -> values
		for(int i = 0; i < this->size(); ++i){ //i is rgs here
			dq[i] = std::vector(*this[i].size());
			for(int j = 0; j < *this[i].size(); ++j){ //j is q's here
				dq[i][j] = std::vector(*this[i][j].size());
				for(k = 0; k < *this[i][j].size()){ //k is dinuc
					int map_q = 0; //maximum a posteriori q
					int best_posterior = 0;
					for(int possible = 0; possible < MAXQ+1; possible++){
						int diff = std::abs(prior[i][j] - possible);
						long double prior_prob = normal_prior[diff];
						long double p = recalibrateutils::q_to_p(possible);
						long double loglike = log_binom_pmf(*this[i][j][k][0], *this[i][j][k][1], p);
						long double posterior = prior_prob + loglike;
						if(posterior > best_posterior){
							map_q = possible;
							best_posterior = posterior;
						}
					}
					dq[i][j][k] = map_q - prior[i][j];
				}
			}
		}
		return dq;
	}

	void CCovariateData::consume_read(const CReadData& read, int minscore = 6){
		rgcov.consume_read(read);
		qcov.consume_read(read);
		cycov.consume_read(read);
		dicov.consume_read(read, minscore);
	}

	dq_t CCovariateData::get_dqs(){
		dq_t dq;
		std::vector<long double> expected_errors(this->qcov.size(),0);
		meanq_t meanq(this->qcov.size(),0);
		for(int rg = 0; rg < this->qcov.size(); ++rg){
			for(int q = 0; q < this->qcov[rg].size(); ++q){
				expected_errors[rg] += q_to_p(q) * this->qcov[rg][q][1];
			}
			meanq[rg] = p_to_q(expected_errors[rg] / this->rgcov[rg][1]);
		}
		dq.meanq = meanq;
		dq.rgdq = this->rgcov.delta_q(meanq);
		prior1_t rgprior(this->qcov.size());
		for(int rg = 0; rg < this->qcov.size(); ++rg){
			qprior[rg].push_back(meanq[rg] + dq.rgdq[rg]);
		}
		dq.qscoredq = this->qcov.delta_q(rgprior);
		prior2_t qprior(this->qcov.size())
		for(int rg = 0; rg < this->qcov.size(); ++rg){
			for(int q = 0; q < this->qcov[rg].size(); ++q){
				qprior.push_back(meanq[rg] + dq.qscoredq[rg][q]);
			}
		}
		dq.cycledq = this->cycov.delta_q(qprior);
		dq.dinucdq = this->dicov.delta_q(qprior);
		return dq;
	}



}