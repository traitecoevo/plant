// -*-c++-*-
#ifndef TREE_ODE_TARGET_H_
#define TREE_ODE_TARGET_H_

#include <vector>

namespace ode {

typedef std::vector<double>::const_iterator iterator_const;
typedef std::vector<double>::iterator       iterator;

class OdeTarget {
public:
  virtual ~OdeTarget() {}

  virtual void derivs(double time, iterator_const y, iterator dydt);
  virtual iterator_const set_ode_values(double time, iterator_const it) = 0;
  virtual iterator       ode_values(iterator it) const = 0;
  virtual iterator       ode_rates(iterator it)  const = 0;
  virtual size_t         ode_size() const = 0;

  std::vector<double> r_derivs(double time, std::vector<double> y);
  void r_set_ode_values(double time, std::vector<double> y);
  std::vector<double> r_ode_values() const;
  std::vector<double> r_ode_rates() const;
};

template <typename TargetIterator>
size_t ode_size(TargetIterator first, TargetIterator last) {
  size_t ret = 0;
  while (first != last) {
    ret += first->ode_size();
    first++;
  }
  return ret;
}

template <typename TargetIterator>
iterator_const set_ode_values(TargetIterator first, TargetIterator last,
			      double time, iterator_const it) {
  while (first != last) {
    it = first->set_ode_values(time, it);
    first++;
  }
  return it;
}

template <typename TargetIterator>
iterator ode_values(TargetIterator first, TargetIterator last,
		    iterator it) {
  while (first != last) {
    it = first->ode_values(it);
    first++;
  }
  return it;
}

template <typename TargetIterator>
iterator ode_rates(TargetIterator first, TargetIterator last,
		   iterator it) {
  while (first != last) {
    it = first->ode_rates(it);
    first++;
  }
  return it;
}

}

#endif
