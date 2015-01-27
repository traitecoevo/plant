#include <tree2.h>

// [[Rcpp::export]]
SEXP oderunner_plant_size(const tree2::ode::Runner<tree2::tools::PlantRunner>& obj) {
  return obj.obj.plant.r_get_vars_size();
}
