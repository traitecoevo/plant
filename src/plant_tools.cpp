#include <plant.h>
#include <plant/uniroot.h>

namespace plant {
namespace tools {

Environment fixed_environment(double canopy_openness,
                              double size_max) {
  std::vector<double> x = {0, size_max/2.0, size_max};
  std::vector<double> y = {canopy_openness, canopy_openness, canopy_openness};
  interpolator::Interpolator env;
  env.init(x, y);
  Parameters<FF16_Strategy> p;
  Environment ret(make_environment(p));
  ret.light_environment = env;
  return ret;
}

}
}

//' Create a light environment where light levels are constant down
//' the canopy.
//'
//' @title Create fixed light environment
//' @param canopy_openness Index of canopy openness (on 0,1)
//' @param size_max Maximum size.  The default (150) should be big
//' enough for most uses.
//' @export
//' @author Rich FitzJohn
// [[Rcpp::export]]
plant::Environment fixed_environment(double canopy_openness,
                                     double size_max=150.0) {
  return plant::tools::fixed_environment(canopy_openness, size_max);
}


// [[Rcpp::export]]
plant::Internals
FF16_oderunner_plant_internals(
  const plant::ode::Runner<plant::tools::PlantRunner<plant::FF16_Strategy>>& obj) {
  return obj.obj.plant.r_internals();
}


// Technical debt: (See RcppR6 #23 and plant #164)
// [[Rcpp::export]]
double FF16_lcp_whole_plant(plant::Plant<plant::FF16_Strategy> p) {
  return plant::tools::lcp_whole_plant(p);
}
