#include "patch.h"

namespace model {

template <>
void Patch<CohortTop>::r_step() {
  r_step_deterministic();
}

SEXP patch(Rcpp::CppClass individual, Parameters p) {
  std::string individual_type =
    util::rcpp_class_demangle(Rcpp::as<std::string>(individual));
  SEXP ret = R_NilValue;
  if (individual_type == "Plant") {
    Patch<Plant> obj(p);
    ret = Rcpp::wrap(obj);
  } else if (individual_type == "CohortDiscrete") {
    Patch<CohortDiscrete> obj(p);
    ret = Rcpp::wrap(obj);
  } else if (individual_type == "CohortTop") {
    Patch<CohortTop> obj(p);
    ret = Rcpp::wrap(obj);
  } else {
    ::Rf_error("Cannot make Patch of %s", individual_type.c_str());
  }
  return ret;
}

}
