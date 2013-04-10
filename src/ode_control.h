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
  void set_eps_abs(double x);
  void set_eps_rel(double x);
  void set_a_y(double x);
  void set_a_dydt(double x);

  double adjust_step_size(size_t dim, unsigned int ord, double step_size,
			  const std::vector<double> &y,
			  const std::vector<double> &yerr,
			  const std::vector<double> &yp);
  double errlevel(double y, double dydt, double step_size);

  bool step_size_shrank() const;
private:
  double eps_abs, eps_rel, a_y, a_dydt;
  bool step_size_shrank_; // TODO: inconsistent name.
};

}

#endif
