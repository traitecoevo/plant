// -*-c++-*-
#ifndef TREE_EBT_H_
#define TREE_EBT_H_

#include "patch.h"
#include "cohort_schedule.h"
// #include "ode_solver.h"

namespace model {

class EBT {
public:
  EBT(Parameters  p);
  EBT(Parameters *p);

  Patch<CohortTop> r_patch() const;

private:
  Patch<CohortTop> patch;
  // ode::Solver<Patch <CohortTop> > ode_solver;
  CohortSchedule schedule;
};

}

#endif
