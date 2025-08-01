// -*-c++-*-
#ifndef PLANT_GET_STATE_H_
#define PLANT_GET_STATE_H_

#include <plant/parameters.h>
#include <plant/node_schedule.h>
#include <plant/scm.h>
#include <plant/stochastic_patch_runner.h>
#include <Rcpp.h>

namespace plant {

  // TODO cleanup scm

// stochastic model:
template <typename T, typename E>
Rcpp::NumericMatrix::iterator get_state(const Individual<T,E>& plant,
                                        Rcpp::NumericMatrix::iterator it) {
  // TODO: this should work (also up in get_state(Node<T,E>, ...)).
  // return plant.ode_state(it);
  std::vector<double> tmp = ode::r_ode_state(plant);
  return std::copy(tmp.begin(), tmp.end(), it);
}

template <typename T, typename E>
Rcpp::NumericMatrix get_state(const StochasticSpecies<T,E>& species) {
  typedef Individual<T,E> individual_type;
  size_t ode_size = individual_type::ode_size(), np = species.size_individuals();
  Rcpp::NumericMatrix ret(static_cast<int>(ode_size), np);

  Rcpp::NumericMatrix::iterator it = ret.begin();
  for (size_t i = 0; i < np; ++i) {
    it = get_state(species.r_individual_at(i), it);
  }
  ret.attr("dimnames") = Rcpp::List::create(individual_type::ode_names(), R_NilValue);
  ret.attr("is_alive") = Rcpp::wrap(species.r_is_alive());
  return ret;
}

template <typename T, typename E>
Rcpp::List get_state(const StochasticPatch<T,E>& patch) {
  Rcpp::List ret;
  for (size_t i = 0; i < patch.size(); ++i) {
    ret.push_back(get_state(patch.at_species(i)));
  }
  return ret;
}

template <typename T, typename E>
Rcpp::List get_state(const StochasticPatchRunner<T,E>& obj) {
  using namespace Rcpp;
  const StochasticPatch<T,E>& patch = obj.r_patch();
  return List::create(_["time"] = obj.time(),
                      _["species"] = get_state(patch));
                      // _["env"] = get_state(patch.r_environment(), obj.time()));
}

}

#endif
