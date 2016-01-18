#include <plant.h>
#include <plant/uniroot.h>

namespace plant {
namespace tools {

Environment fixed_environment(double canopy_openness,
                              double height_max) {
  std::vector<double> x = {0, height_max/2.0, height_max};
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

// TOOD: Now that this returns all the variables, it probably needs a
// change of name.  However, it's only used internally so it's not
// that big of a deal.
// [[Rcpp::export]]
plant::PlantPlus_internals
oderunner_plant_size(const plant::ode::Runner<plant::tools::PlantRunner>& obj) {
  return obj.obj.plant.r_internals();
}

//' Create a light environment where light levels are constant down
//' the canopy.
//'
//' @title Create fixed light environment
//' @param canopy_openness Index of canopy openness (on 0,1)
//' @param height_max Maximum height.  The default (150) should be big
//' enough for most uses.
//' @export
//' @author Rich FitzJohn
// [[Rcpp::export]]
plant::Environment fixed_environment(double canopy_openness,
                                     double height_max=150.0) {
  return plant::tools::fixed_environment(canopy_openness, height_max);
}


// Technical debt: (See RcppR6 #23 and plant #164)

// [[Rcpp::export]]
double FF16_lcp_whole_plant(plant::PlantPlus<plant::FF16_Strategy> p) {
  return plant::tools::lcp_whole_plant(p);
}
// [[Rcpp::export]]
double FF16r_lcp_whole_plant(plant::PlantPlus<plant::FF16r_Strategy> p) {
  return plant::tools::lcp_whole_plant(p);
}
