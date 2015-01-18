// -*-c++-*-
#ifndef TREE_ODE_INTERFACE_H_
#define TREE_ODE_INTERFACE_H_

namespace ode {

// These are utilities designed to make it more pleasant to work with
// ode objects.  The first in-place functions work with the
// boost::odeint interface, and the types are just to reuce typing.
typedef std::vector<double>        state_type;
typedef state_type::const_iterator const_iterator;
typedef state_type::iterator       iterator;

template <typename T>
size_t ode_size(const T& obj) {
  return obj.ode_size();
}

// For models that do not care about time, we leave time out here.
template <typename T>
void set_ode_values(T& obj, const state_type& y) {
  util::check_length(y.size(), obj.ode_size());
  obj.set_ode_values(y.begin());
}

template <typename T>
void set_ode_values(T& obj, const state_type& y, double time) {
  util::check_length(y.size(), obj.ode_size());
  obj.set_ode_values(y.begin(), time);
}

// NOTE: Might need to also have a non-const version here?
template <typename T>
void ode_values(const T& obj, state_type& values) {
  obj.ode_values(values.begin());
}

template <typename T>
void ode_rates(const T& obj, state_type& dydt) {
  obj.ode_rates(dydt.begin());
}

// These out-of-place versions are useful for interfacing with R.
template <typename T>
state_type r_ode_values(const T& obj) {
  state_type values(obj.ode_size());
  ode_values(obj, values);
  return values;
}

template <typename T>
state_type r_ode_rates(const T& obj) {
  state_type dydt(obj.ode_size());
  ode_rates(obj, dydt);
  return dydt;
}

// This is the type expected by odeint, but with the first argument
// being a reference to the object.  It's going to need to be a
// mutable reference in our case because the derivative function is
// going to need to set and change a number of things within the
// object.  All that needs doing after this is binding with a lambda
// function:
//
//   auto f = [this] (const ode::state_type& y, state_type& dydt,
//                    const double time) -> void {
//     derivs(this, y, dydt, time);
//   };
template <typename T>
void derivs(T& obj, const state_type& y, state_type& dydt,
            const double time) {
  set_ode_values(obj, y, time);
  ode_rates(obj, dydt);
}

}

#endif
