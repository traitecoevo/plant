#include "cohort_discrete.h"

namespace model {

CohortDiscrete::CohortDiscrete(Strategy s)
  : Plant(s),
    n_individuals(1) {
}

CohortDiscrete::CohortDiscrete(Strategy *s)
  : Plant(s),
    n_individuals(1) {
}

CohortDiscrete::CohortDiscrete(Strategy s, int n_individuals_)
  : Plant(s),
    n_individuals(n_individuals_) {
  if (n_individuals < 1)
    Rcpp::stop("Cannot create a cohort with less than 1 individual");
}

CohortDiscrete::CohortDiscrete(Strategy *s, int n_individuals_)
  : Plant(s),
    n_individuals(n_individuals_) {
  if (n_individuals < 1)
    Rcpp::stop("Cannot create a cohort with less than 1 individual");
}

double CohortDiscrete::leaf_area() const {
  return n_individuals * Plant::leaf_area();
}

double CohortDiscrete::leaf_area_above(double z) const {
  return n_individuals * Plant::leaf_area_above(z);
}

int CohortDiscrete::offspring() {
  return n_individuals * Plant::offspring();
}

bool CohortDiscrete::died() {
  const double p_died = mortality_probability();
  if (n_individuals > 1)
    n_individuals -= static_cast<int>(Rf_rbinom(n_individuals, p_died));
  else if (n_individuals == 1) // ensures same RNG behaviour as Plant
    n_individuals = unif_rand() < p_died ? 0 : 1;
  set_mortality(0.0);
  if (n_individuals < 0)
    Rcpp::stop("Somehow we have negative individuals!");
  return n_individuals == 0;
}

int CohortDiscrete::get_n_individuals() const {
  return n_individuals;
}

void CohortDiscrete::set_n_individuals(int n) {
  if (n < 1)
    Rcpp::stop("n must be positive");
  n_individuals = n;
}

size_t CohortDiscrete::state_size() const {
  return ode_size() + 1;
}

CohortDiscrete::state::iterator
CohortDiscrete::get_state(CohortDiscrete::state::iterator it) const {
  it = Plant::get_state(it);
  *it++ = static_cast<double>(n_individuals);
  return it;
}

CohortDiscrete::state::const_iterator
CohortDiscrete::set_state(CohortDiscrete::state::const_iterator it) {
  it = Plant::set_state(it);
  n_individuals = static_cast<int>(*it++);
  return it;
}

}
