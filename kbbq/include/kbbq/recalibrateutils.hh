#ifndef RECALIBRATEUTILS_H
#define RECALIBRATEUTILS_H

#include <vector>
#include <string>
#include "readutils.hh"
#include "covariateutils.hh"
#include "bloom.hh"

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

bloom::bloomary_t subsample_kmers(std::string filename, int ksize, double alpha);

}

#endif