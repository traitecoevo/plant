// -*-c++-*-
#ifndef TREE_ODE_TARGET_H_
#define TREE_ODE_TARGET_H_

#include <vector>

namespace ode {

typedef std::vector<double>::const_iterator iter_const;
typedef std::vector<double>::iterator       iter;

class OdeTarget {
public:
  virtual ~OdeTarget() {}

  // virtual void derivs(double time, iter_const y, iter dydt) = 0;
  virtual iter_const ode_values_set(iter_const it, bool &changed) = 0;
  virtual iter       ode_values(iter it) const = 0;
  virtual iter       ode_rates(iter it)  const = 0;
  virtual size_t     ode_size() const = 0;

  bool r_ode_values_set(std::vector<double> y);
  std::vector<double> r_ode_values() const;
  std::vector<double> r_ode_rates() const;
};

}


#endif
