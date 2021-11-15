// -*-c++-*-
#ifndef PLANT_GET_AUX_H_
#define PLANT_GET_AUX_H_

#include <plant/parameters.h>
#include <plant/cohort_schedule.h>
#include <plant/scm.h>
#include <Rcpp.h>

namespace plant {

template <typename T, typename E>
Rcpp::NumericMatrix::iterator get_aux(const Cohort<T,E>& cohort,
                                        Rcpp::NumericMatrix::iterator it) {
  std::vector<double> tmp = ode::r_ode_aux(cohort);
  return std::copy(tmp.begin(), tmp.end(), it);
}

template <typename T, typename E>
Rcpp::NumericMatrix get_aux(const Species<T,E>& species) {
  typedef Cohort<T,E> cohort_type;
  size_t aux_size = species.strategy_aux_size(), np = species.size();
  Rcpp::NumericMatrix ret(static_cast<int>(aux_size), np + 1); // +1 is seed
  Rcpp::NumericMatrix::iterator it = ret.begin();
  for (size_t i = 0; i < np; ++i) {
    it = get_aux(species.r_cohort_at(i), it);
  }
  it = get_aux(species.r_new_cohort(), it);
  ret.attr("dimnames") =
    Rcpp::List::create(species.aux_names(), R_NilValue);
  return ret;
}

template <typename T, typename E>
Rcpp::List get_aux(const Patch<T,E>& patch) {
  Rcpp::List ret;
  for (size_t i = 0; i < patch.size(); ++i) {
    ret.push_back(get_aux(patch.at(i)));
  }
  return ret;
}

template <typename T, typename E>
Rcpp::List get_aux(const SCM<T,E>& scm) {
  using namespace Rcpp;
  const Patch<T,E>& patch = scm.r_patch();
  return List::create(_["species"] = get_aux(patch));
}

}

#endif