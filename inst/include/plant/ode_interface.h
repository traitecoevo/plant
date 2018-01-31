// -*-c++-*-
#ifndef PLANT_PLANT_ODE_INTERFACE_H_
#define PLANT_PLANT_ODE_INTERFACE_H_

namespace plant {
namespace ode {

// These are utilities designed to make it more pleasant to work with
// ode objects.  The first in-place functions work with the
// boost::odeint interface, and the types are just to reduce typing.
typedef std::vector<double>        state_type;
typedef state_type::const_iterator const_iterator;
typedef state_type::iterator       iterator;

// By default, we assume that systems are time homogeneous; systems
// that provide an `ode_time` function will be treated differently.
template <typename T>
class needs_time {
  typedef char true_type;
  typedef long false_type;
  template <typename C> static true_type test(decltype(&C::ode_time)) ;
  template <typename C> static false_type test(...);
public:
  enum { value = sizeof(test<T>(0)) == sizeof(true_type) };
};

// The recursive interface
template <typename ForwardIterator>
size_t ode_size(ForwardIterator first, ForwardIterator last) {
  size_t ret = 0;
  while (first != last) {
    ret += first->ode_size();
    ++first;
  }
  return ret;
}

template <typename ForwardIterator>
const_iterator set_ode_state(ForwardIterator first, ForwardIterator last,
                             const_iterator it) {
  while (first != last) {
    it = first->set_ode_state(it);
    ++first;
  }
  return it;
}

template <typename ForwardIterator>
iterator ode_state(ForwardIterator first, ForwardIterator last,
                   iterator it) {
  while (first != last) {
    it = first->ode_state(it);
    ++first;
  }
  return it;
}

template <typename ForwardIterator>
iterator ode_rates(ForwardIterator first, ForwardIterator last,
                   iterator it) {
  while (first != last) {
    it = first->ode_rates(it);
    ++first;
  }
  return it;
}

template <typename T>
typename std::enable_if<needs_time<T>::value, double>::type
ode_time(const T& obj) {
  return obj.ode_time();
}

template <typename T>
typename std::enable_if<!needs_time<T>::value, double>::type
ode_time(const T& /* obj */) {
  return 0.0;
}

namespace internal {
template <typename T>
typename std::enable_if<needs_time<T>::value, void>::type
set_ode_state(T& obj, const state_type& y, double time) {
  obj.set_ode_state(y.begin(), time);
}

template <typename T>
typename std::enable_if<!needs_time<T>::value, void>::type
set_ode_state(T& obj, const state_type& y, double /* time */) {
  obj.set_ode_state(y.begin());
}
}

template <typename T>
void derivs(T& obj, const state_type& y, state_type& dydt,
            const double time) {
  internal::set_ode_state(obj, y, time);
  obj.ode_rates(dydt.begin());
}

template <typename T>
state_type r_derivs(T& obj, const state_type& y, const double time) {
  state_type dydt(obj.ode_size());
  derivs(obj, y, dydt, time);
  return dydt;
}

// These out-of-place versions are useful for interfacing with R.
template <typename T>
typename std::enable_if<needs_time<T>::value, void>::type
r_set_ode_state(T& obj, const state_type& y, double time) {
  util::check_length(y.size(), obj.ode_size());
  obj.set_ode_state(y.begin(), time);
}

template <typename T>
typename std::enable_if<!needs_time<T>::value, void>::type
r_set_ode_state(T& obj, const state_type& y) {
  util::check_length(y.size(), obj.ode_size());
  obj.set_ode_state(y.begin());
}

template <typename T>
typename std::enable_if<needs_time<T>::value, double>::type
r_ode_time(const T& obj) {
  return obj.ode_time();
}

template <typename T>
typename std::enable_if<!needs_time<T>::value, double>::type
r_ode_time(const T& /* obj */) {
  return 0.0;
}

template <typename T>
state_type r_ode_state(const T& obj) {
  state_type values(obj.ode_size());
  obj.ode_state(values.begin());
  return values;
}

template <typename T>
state_type r_ode_rates(const T& obj) {
  state_type dydt(obj.ode_size());
  obj.ode_rates(dydt.begin());
  return dydt;
}

}
}

#endif
