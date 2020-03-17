#ifndef RECALIBRATEUTILS_H
#define RECALIBRATEUTILS_H

#include <vector>
#include "readutils.h"

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


}

#endif