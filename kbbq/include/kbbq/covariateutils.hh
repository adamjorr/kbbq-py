#ifndef COVARIATEUTILS_H
#define COVARIATEUTILS_H
#define MAXQ 93

#include <vector>
#include <utility>
#include <cmath>
#include <random>
#include "readutils.hh"
#include "recalibrateutils.hh"

namespace covariateutils{

typedef meanq_t std::vector<int>;
typedef rgdq_t std::vector<int>;
typedef qscoredq_t std::vector<std::vector<int>>;
typedef cycledq_t std::vector<std::vector<std::vector<int>>>;
typedef dinucdq_t std::vector<std::vector<std::vector<int>>>;

struct dq_t
{
	meanq_t meanq;
	globaldq_t globaldq;
	rgdq_t qscoredq;
	cycledq_t cycledq;
	dinucdq_t dinucdq;
};

std::vector<long double> normal_prior(MAXQ+1);

for(int i = 0; i < MAXQ+1; ++i){
	//we might need to do something here to handle over/underflow but YOLO
	normal_prior[i] = std::log(.9l * std::expl(-(((long double)i/.5)**2)/2));
}

typedef covariate_t std::array<int, 2>;

class CCovariate: public std::vector<covariate_t>
{
public:
	CCovariate(size_t len): std::vector<covariate_t>(len) {}
	void increment(size_t idx, covariate_t value);
};

class CRGCovariate: public CCovariate
{
public:
	CRGCovariate(size_t len): CCovariate(len){}
	void consume_read(const CReadData& read);
	rgdq_t delta_q(int prior, int max);
}

class CQCovariate: public std::vector<CCovariate>
{
public:
	CQCovariate(size_t rgs, size_t qlen): std::vector<CCovariate>(rgs, CCovariate(qlen)){}
	void consume_read(const CReadData& read);
	dinucdq_t delta_q(int prior, int max);
}

typedef cycle_t std::array<CCovariate,2>;

//The first Covariate is for fwd reads, the 2nd Covariate is for reverse reads.
class CCycleCovariate: public std::vector<std::vector<cycle_t>>
{
public:
	CCycleCovariate(size_t rgs, size_t qlen, size_t cylen):
		std::vector<std::vector<cycle_t>>(rgs, std::vector<cycle_t>(qlen, cycle_t(CCovariate(cylen),CCovariate(cylen))))
		{}
	void consume_read(const CReadData& read);
	cycledq_t delta_q(int prior, int max);
}

class CDinucCovariate: public std::vector<std::vector<CCovariate>>
{
public:
	CDinucCovariate(size_t rgs, size_t qlen, size_t dilen):
		std::vector<std::vector<CCovariate>>(rgs, std::vector<CCovariate>(qlen, CCovariate(dilen)))
		{}
	void consume_read(const CReadData& read, int minscore = 6);
	dinucdq_t delta_q(int prior, int max);
}

class CCovariateData():
{
protected:
	CRGCovariate rgcov;
	CQCovariate qcov;
	CCycleCovariate cycov;
	CDinucCovariate dicov;
public:
	CCovariateData(){};
	void consume_read(const CReadData& read, int minscore = 6);
	dq_t get_dqs(covariateutils::CCovariateData);
}

}
#endif