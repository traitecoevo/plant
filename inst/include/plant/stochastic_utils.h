// -*-c++-*-
#ifndef PLANT_PLANT_STOCHASTIC_UTILS_H_
#define PLANT_PLANT_STOCHASTIC_UTILS_H_

#include <plant/node_schedule.h>
#include <plant/scm_utils.h>

namespace plant {

template <typename Parameters>
NodeSchedule make_empty_stochastic_schedule(const Parameters& p) {
  NodeSchedule ret(p.size());
  ret.r_set_max_time(p.max_patch_lifetime);
  return ret;
}

}

#endif
