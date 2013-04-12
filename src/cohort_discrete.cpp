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

CohortDiscrete::CohortDiscrete(Strategy s, int n_individuals)
  : Plant(s),
    n_individuals(n_individuals) {
  if ( n_individuals < 1 ) // TODO: replace with throw
    ::Rf_error("Cannot create a cohort with less than 1 individual");
}

CohortDiscrete::CohortDiscrete(Strategy *s, int n_individuals)
  : Plant(s),
    n_individuals(n_individuals) {
  if ( n_individuals < 1 ) // TODO: replace with throw
    ::Rf_error("Cannot create a cohort with less than 1 individual");
}

double CohortDiscrete::leaf_area_above(double z) const {
  return n_individuals * Plant::leaf_area_above(z);
}

int CohortDiscrete::offspring() {
  return n_individuals * Plant::offspring();
}

bool CohortDiscrete::died() {
  if ( n_individuals > 1 )
    n_individuals -= Rf_rbinom(n_individuals, mortality());
  else if ( n_individuals == 1 ) // ensures same RNG behaviour as Plant
    n_individuals = unif_rand() < mortality() ? 0 : 1;
  mortality_reset();
  if ( n_individuals < 0 )
    Rf_error("Somehow we have negative individuals!");
  return n_individuals == 0;
}

int CohortDiscrete::r_n_individuals() const {
  return n_individuals;
}

void CohortDiscrete::r_set_n_individuals(int n) {
  n_individuals = n;
}

}
