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
}

CohortDiscrete::CohortDiscrete(Strategy *s, int n_individuals)
  : Plant(s),
    n_individuals(n_individuals) {
}

double CohortDiscrete::leaf_area_above(double z) const {
  return n_individuals * leaf_area_above(z);
}

int CohortDiscrete::offspring() {
  return n_individuals * offspring();
}

bool CohortDiscrete::died() {
  n_individuals -= Rf_rbinom(n_individuals, mortality());
  mortality_reset();
  if ( n_individuals < 0 )
    Rf_error("Somehow we have negative individuals!");
  return n_individuals == 0;
}

}
