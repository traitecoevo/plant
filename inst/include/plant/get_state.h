// -*-c++-*-
#ifndef PLANT_GET_STATE_H_
#define PLANT_GET_STATE_H_

#include <plant/parameters.h>
#include <plant/cohort_schedule.h>
#include <plant/ebt.h>
#include <plant/stochastic_patch_runner.h>
#include <Rcpp.h>

namespace plant {

// TODO: Might be worth moving this into the main interface?  Not sure
// about that though as it revolves pretty closely around Rcpp types
// and is kind of separate to the rest of the model.  It won't matter
// though; if it does move in then we just adjust the yml.
template <typename T>
Rcpp::NumericMatrix::iterator get_state(const Cohort<T>& cohort,
                                        Rcpp::NumericMatrix::iterator it) {
  std::vector<double> tmp = ode::r_ode_state(cohort);
  return std::copy(tmp.begin(), tmp.end(), it);
}

template <typename T>
Rcpp::NumericMatrix get_state(const Species<T>& species) {
  typedef Cohort<T> cohort_type;
  const size_t ode_size = cohort_type::ode_size(), np = species.size();
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

template <typename T>
Rcpp::List get_state(const Patch<T>& patch) {
  Rcpp::List ret;
  for (size_t i = 0; i < patch.size(); ++i) {
    ret.push_back(get_state(patch.at(i)));
  }
  return ret;
}

inline Rcpp::NumericMatrix get_state(const Environment environment) {
  using namespace Rcpp;
  NumericMatrix xy = environment.light_environment.r_get_xy();
  Rcpp::CharacterVector colnames =
    Rcpp::CharacterVector::create("height", "canopy_openness");
  xy.attr("dimnames") = Rcpp::List::create(R_NilValue, colnames);
  return xy;
}

template <typename T>
Rcpp::List get_state(const EBT<T>& ebt) {
  using namespace Rcpp;
  const Patch<T>& patch = ebt.r_patch();
  return List::create(_["time"] = ebt.time(),
                      _["species"] = get_state(patch),
                      _["light_env"] = get_state(patch.r_environment()));
}

// stochastic model:
template <typename T>
Rcpp::NumericMatrix::iterator get_state(const Plant<T>& plant,
                                        Rcpp::NumericMatrix::iterator it) {
  // TODO: this should work (also up in get_state(Cohort<T>, ...)).
  // return plant.ode_state(it);
  std::vector<double> tmp = ode::r_ode_state(plant);
  return std::copy(tmp.begin(), tmp.end(), it);
}

template <typename T>
Rcpp::NumericMatrix get_state(const StochasticSpecies<T>& species) {
  typedef Plant<T> plant_type;
  const size_t ode_size = plant_type::ode_size(), np = species.size_plants();
  Rcpp::NumericMatrix ret(static_cast<int>(ode_size), np);

  Rcpp::NumericMatrix::iterator it = ret.begin();
  for (size_t i = 0; i < np; ++i) {
    it = get_state(species.r_plant_at(i), it);
  }
  ret.attr("dimnames") =
    Rcpp::List::create(plant_type::ode_names(), R_NilValue);
  ret.attr("is_alive") = Rcpp::wrap(species.r_is_alive());
  return ret;
}

template <typename T>
Rcpp::List get_state(const StochasticPatch<T>& patch) {
  Rcpp::List ret;
  for (size_t i = 0; i < patch.size(); ++i) {
    ret.push_back(get_state(patch.at(i)));
  }
  return ret;
}

template <typename T>
Rcpp::List get_state(const StochasticPatchRunner<T>& obj) {
  using namespace Rcpp;
  const StochasticPatch<T>& patch = obj.r_patch();
  return List::create(_["time"] = obj.time(),
                      _["species"] = get_state(patch),
                      _["light_env"] = get_state(patch.r_environment()));
};

}

#endif
