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

  void step();
  void run_next();

  double get_time() const;

  Patch<CohortTop> r_patch() const;
  CohortSchedule r_cohort_schedule() const;
  void r_set_cohort_schedule(CohortSchedule x);

private:
  Patch<CohortTop> patch;
  ode::Solver<Patch <CohortTop> > ode_solver;
  CohortSchedule schedule;
};

}

#endif
