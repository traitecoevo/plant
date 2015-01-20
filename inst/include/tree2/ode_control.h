// -*-c++-*-
#ifndef TREE_ODE_CONTROL_H_
#define TREE_ODE_CONTROL_H_

#include <vector>
#include <cstddef>

namespace ode {

struct OdeControl {
  typedef std::vector<double> state_type;
  OdeControl();
  OdeControl(double tol_abs_, double tol_rel_,
	     double a_y_, double a_dydt_,
	     double step_size_min_, double step_size_max_,
	     double step_size_initial_);

  double adjust_step_size(size_t dim, size_t ord, double step_size,
			  const state_type& y,
			  const state_type& yerr,
			  const state_type& yp);
  double errlevel(double y, double dydt, double step_size) const;
  bool step_size_shrank() const;

  double tol_abs, tol_rel, a_y, a_dydt;
  double step_size_min, step_size_max, step_size_initial;
  bool last_step_size_shrank;
};

}

#endif
