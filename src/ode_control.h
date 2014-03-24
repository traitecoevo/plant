// -*-c++-*-
#ifndef TREE_ODE_CONTROL_H_
#define TREE_ODE_CONTROL_H_

#include <Rcpp.h>
#include <vector>
#include <cstddef>
#include "lookup.h"

// There is also a "scaled" control type which has a vector of scaling
// parameters along the model parameters; this could be added very
// easily here with an extra slot in the class.

namespace ode {

class OdeControl : public util::Lookup {
public:
  OdeControl();
  OdeControl(double tol_abs_, double tol_rel_,
	     double a_y_, double a_dydt_,
	     double step_size_min_, double step_size_max_);

  double adjust_step_size(size_t dim, size_t ord, double step_size,
			  const std::vector<double> &y,
			  const std::vector<double> &yerr,
			  const std::vector<double> &yp);
  double errlevel(double y, double dydt, double step_size) const;

  bool step_size_shrank() const;
private:
  void do_build_lookup();
  double tol_abs, tol_rel, a_y, a_dydt;
  double step_size_min, step_size_max;
  bool last_step_size_shrank;
};

}

RCPP_EXPOSED_CLASS_NODECL(ode::OdeControl)

#endif
