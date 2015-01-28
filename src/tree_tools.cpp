#include <tree2.h>
#include <tree2/uniroot.h>

namespace tree2 {
namespace tools {

Environment fixed_environment(double canopy_openness,
			      double height_max) {
  std::vector<double> x = {0, height_max/2.0, height_max};
  std::vector<double> y = {canopy_openness, canopy_openness, canopy_openness};
  interpolator::Interpolator env;
  env.init(x, y);
  Parameters p;
  Environment ret(p);
  ret.light_environment = env;
  return ret;
}

double lcp_whole_plant(Plant p) {
  auto target = [&] (double x) mutable -> double {
    Environment env = fixed_environment(x);
    p.compute_vars_phys(env);
    return p.net_production();
  };

  const double f1 = target(1.0);
  if (f1 < 0.0) {
    return NA_REAL;
  } else {
    const double tol = p.control().plant_seed_tol;
    const size_t max_iterations = p.control().plant_seed_iterations;
    return util::uniroot(target, 0.0, 1.0, tol, max_iterations);
  }
}

}
}

// [[Rcpp::export]]
SEXP oderunner_plant_size(const tree2::ode::Runner<tree2::tools::PlantRunner>& obj) {
  return obj.obj.plant.r_get_vars_size();
}

// [[Rcpp::export]]
tree2::Environment fixed_environment(double canopy_openness,
				     double height_max=150.0) {
  return tree2::tools::fixed_environment(canopy_openness, height_max);
}


// [[Rcpp::export]]
double lcp_whole_plant(tree2::Plant p) {
  return tree2::tools::lcp_whole_plant(p);
}

// // [--[Rcpp::export(lcp_whole_plant)]]
// double r_lcp_whole_plant(tree2::Plant p) {
//   return tree2::tools::lcp_whole_plant(p);
// }
