#include <tree.h>
#include <tree/uniroot.h>

namespace tree {
namespace tools {

Environment fixed_environment(double canopy_openness,
			      double height_max) {
  std::vector<double> x = {0, height_max/2.0, height_max};
  std::vector<double> y = {canopy_openness, canopy_openness, canopy_openness};
  interpolator::Interpolator env;
  env.init(x, y);
  FFW16_Parameters p;
  Environment ret(p);
  ret.light_environment = env;
  return ret;
}

double lcp_whole_plant(FFW16_PlantPlus p) {
  auto target = [&] (double x) mutable -> double {
    Environment env = fixed_environment(x);
    p.compute_vars_phys(env);
    return p.net_mass_production_dt();
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

// TOOD: Now that this returns all the variables, it probably needs a
// change of name.  However, it's only used internally so it's not
// that big of a deal.
// [[Rcpp::export]]
tree::FFW16_PlantPlus_internals
oderunner_plant_size(const tree::ode::Runner<tree::tools::PlantRunner>& obj) {
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
tree::Environment fixed_environment(double canopy_openness,
				     double height_max=150.0) {
  return tree::tools::fixed_environment(canopy_openness, height_max);
}

//' Compute the whole plant light compensation point for a single
//' plant.
//' @title Whole plant light compensation point
//' @param p A \code{Plant}, with strategy, height, etc set.
//' @export
//' @author Rich FitzJohn
// [[Rcpp::export]]
double lcp_whole_plant(tree::FFW16_PlantPlus p) {
  return tree::tools::lcp_whole_plant(p);
}
