// -*-c++-*-
#ifndef TREE_CONTROL_H_
#define TREE_CONTROL_H_

#include <Rcpp.h>

#include "lookup.h"
#include "util.h"

namespace model {

class Control : public util::Lookup {
public:
  typedef util::PtrWrapper<Control> ptr;
  Control();
  Control(Rcpp::List x);

  bool   plant_assimilation_over_distribution;
  double plant_assimilation_tol;
  int    plant_assimilation_iterations;

  double plant_seed_tol;
  int    plant_seed_iterations;

  double cohort_gradient_eps;
  bool   cohort_gradient_richardson;
  int    cohort_gradient_richardson_depth;

private:
  double _plant_assimilation_over_distribution;
  double _plant_assimilation_iterations;

  double _plant_seed_iterations;

  double _cohort_gradient_richardson;
  double _cohort_gradient_richardson_depth;

  void do_build_lookup();
  void reset();
  void set_parameters_post_hook();
};

}

RCPP_EXPOSED_CLASS(model::Control)

#endif
