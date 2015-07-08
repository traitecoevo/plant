// -*-c++-*-
#ifndef PLANT_PLANT_UTILS_H_
#define PLANT_PLANT_UTILS_H_

#include <plant/parameters.h>
#include <plant/cohort_schedule.h>
#include <plant/ebt.h>
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
  Rcpp::NumericMatrix ret(static_cast<int>(species.r_seed().ode_size()),
			  species.size() + 1); // +1 is seed
  Rcpp::NumericMatrix::iterator it = ret.begin();
  for (size_t i = 0; i < species.size(); ++i) {
    it = get_state(species.r_cohort_at(i), it);
  }
  it = get_state(species.r_seed(), it);

  // Add dimension names to this:
  Rcpp::CharacterVector rownames =
    Rcpp::CharacterVector::create("height", "log_mortality",
                                  "seeds", "log_density");
  // Can probably do this via static_assert?
  if (ret.nrow() != rownames.size()) {
    Rcpp::stop("Hmm - looks like the model has changed...");
  }
  ret.attr("dimnames") = Rcpp::List::create(rownames, R_NilValue);

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

}


#endif
