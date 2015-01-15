// -*-c++-*-
#ifndef TREE_COHORT_H_
#define TREE_COHORT_H_

#include <tree2/environment.h>

namespace tree2 {

template <typename T>
class Cohort {
public:
  typedef typename T::strategy_type strategy_type;
  Cohort(strategy_type s);
  T plant;
};

template <typename T>
Cohort<T>::Cohort(strategy_type s)
  : plant(s) {
}

template <typename T>
Cohort<T> make_cohort(typename Cohort<T>::strategy_type::element_type s) {
  return Cohort<T>(make_strategy_ptr(s));
}

}

#endif
