// -*-c++-*-
#ifndef TREE2_TREE_UTILS_H_
#define TREE2_TREE_UTILS_H_

#include <tree2/parameters.h>
#include <tree2/cohort_schedule.h>

namespace tree2 {

std::vector<double> cohort_schedule_default_times(double max_time);
double cohort_schedule_default_max_time(const Parameters& p);
CohortSchedule cohort_schedule_default(const Parameters& p);

}


#endif
