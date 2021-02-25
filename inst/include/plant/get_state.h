// -*-c++-*-
#ifndef PLANT_GET_STATE_H_
#define PLANT_GET_STATE_H_

#include <plant/parameters.h>
#include <plant/cohort_schedule.h>
#include <plant/scm.h>
#include <plant/stochastic_patch_runner.h>
#include <Rcpp.h>

namespace plant {

// TODO: Might be worth moving this into the main interface?  Not sure
// about that though as it revolves pretty closely around Rcpp types
// and is kind of separate to the rest of the model.  It won't matter
// though; if it does move in then we just adjust the yml.
template <typename T, typename E>
Rcpp::NumericMatrix::iterator get_state(const Cohort<T,E>& cohort,
                                        Rcpp::NumericMatrix::iterator it) {
  std::vector<double> tmp = ode::r_ode_state(cohort);
  return std::copy(tmp.begin(), tmp.end(), it);
}

template <typename T, typename E>
Rcpp::NumericMatrix get_state(const Species<T,E>& species) {
  typedef Cohort<T,E> cohort_type;
  size_t ode_size = cohort_type::ode_size(), np = species.size();
  Rcpp::NumericMatrix ret(static_cast<int>(ode_size), np + 1); // +1 is seed
  Rcpp::NumericMatrix::iterator it = ret.begin();
  for (size_t i = 0; i < np; ++i) {
    it = get_state(species.r_cohort_at(i), it);
  }
  it = get_state(species.r_seed(), it);
  ret.attr("dimnames") =
    Rcpp::List::create(cohort_type::ode_names(), R_NilValue);
  return ret;
}

template <typename T, typename E>
Rcpp::List get_state(const Patch<T,E>& patch) {
  Rcpp::List ret;
  for (size_t i = 0; i < patch.size(); ++i) {
    ret.push_back(get_state(patch.at(i)));
  }
  return ret;
}

inline Rcpp::NumericMatrix get_state(const Environment environment) {
  // Empty vector
  std::vector<std::vector<double>> xy;
  return Rcpp::wrap(util::to_rcpp_matrix(xy));
}

template <typename T, typename E>
Rcpp::List get_state(const SCM<T,E>& scm) {
  using namespace Rcpp;
  const Patch<T,E>& patch = scm.r_patch();
  return List::create(_["time"] = scm.time(),
                      _["species"] = get_state(patch),
                      _["env"] = get_state(patch.r_environment()));
}

// stochastic model:
template <typename T, typename E>
Rcpp::NumericMatrix::iterator get_state(const Individual<T,E>& plant,
                                        Rcpp::NumericMatrix::iterator it) {
  // TODO: this should work (also up in get_state(Cohort<T,E>, ...)).
  // return plant.ode_state(it);
  std::vector<double> tmp = ode::r_ode_state(plant);
  return std::copy(tmp.begin(), tmp.end(), it);
}

template <typename T, typename E>
Rcpp::NumericMatrix get_state(const StochasticSpecies<T,E>& species) {
  typedef Individual<T,E> individual_type;
  size_t ode_size = individual_type::ode_size(), np = species.size_plants();
  Rcpp::NumericMatrix ret(static_cast<int>(ode_size), np);

  Rcpp::NumericMatrix::iterator it = ret.begin();
  for (size_t i = 0; i < np; ++i) {
    it = get_state(species.r_plant_at(i), it);
  }
  ret.attr("dimnames") =
    Rcpp::List::create(individual_type::ode_names(), R_NilValue);
  ret.attr("is_alive") = Rcpp::wrap(species.r_is_alive());
  return ret;
}

template <typename T, typename E>
Rcpp::List get_state(const StochasticPatch<T,E>& patch) {
  Rcpp::List ret;
  for (size_t i = 0; i < patch.size(); ++i) {
    ret.push_back(get_state(patch.at(i)));
  }
  return ret;
}

template <typename T, typename E>
Rcpp::List get_state(const StochasticPatchRunner<T,E>& obj) {
  using namespace Rcpp;
  const StochasticPatch<T,E>& patch = obj.r_patch();
  return List::create(_["time"] = obj.time(),
                      _["species"] = get_state(patch),
                      _["env"] = get_state(patch.r_environment()));
};

}

#endif
