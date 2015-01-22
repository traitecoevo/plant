// -*-c++-*-
#ifndef TREE_EBT_UTILS_H_
#define TREE_EBT_UTILS_H_

#include <tree2/cohort_schedule.h>
#include <tree2/parameters.h>

namespace tree2 {

std::vector<double> cohort_schedule_times_default(double max_time);
double cohort_schedule_max_time_default(const Parameters& p);
CohortSchedule cohort_schedule_default(const Parameters& p);
CohortSchedule make_cohort_schedule(const Parameters& p);

}

#endif
