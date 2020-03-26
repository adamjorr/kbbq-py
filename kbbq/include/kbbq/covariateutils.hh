#ifndef COVARIATEUTILS_H
#define COVARIATEUTILS_H

#include <vector>
#include <utility>
#include "readutils.hh"

namespace covariateutils{

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
}

class CQCovariate: public std::vector<CCovariate>
{
public:
	CQCovariate(size_t rgs, size_t qlen): std::vector<CCovariate>(rgs, CCovariate(qlen)){}
	void consume_read(const CReadData& read);
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
}

class CDinucCovariate: public std::vector<std::vector<CCovariate>>
{
public:
	CDinucCovariate(size_t rgs, size_t qlen, size_t dilen):
		std::vector<std::vector<CCovariate>>(rgs, std::vector<CCovariate>(qlen, CCovariate(dilen)))
		{}
	void consume_read(const CReadData& read, int minscore = 6);
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
	void consume_read(const CReadData& read, int minscore);
}

}
#endif