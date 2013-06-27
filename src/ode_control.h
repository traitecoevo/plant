// -*-c++-*-
#ifndef TREE_ODE_CONTROL_H_
#define TREE_ODE_CONTROL_H_

#include <vector>

// There is also a "scaled" control type which has a vector of scaling
// parameters along the model parameters; this could be added very
// easily here with an extra slot in the class.

namespace ode {

class OdeControl {
public:
  OdeControl();
  OdeControl(double tol_abs, double tol_rel,
	     double a_y, double a_dydt,
	     double step_size_min, double step_size_max);

  double adjust_step_size(size_t dim, unsigned int ord, double step_size,
			  const std::vector<double> &y,
			  const std::vector<double> &yerr,
			  const std::vector<double> &yp);
  double errlevel(double y, double dydt, double step_size) const;

  bool step_size_shrank() const;
private:
  double tol_abs, tol_rel, a_y, a_dydt;
  double step_size_min, step_size_max;
  bool last_step_size_shrank;
};

}

#endif
