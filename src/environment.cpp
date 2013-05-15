#include "environment.h"

namespace model {

Environment::Environment(spline::Spline *light_environment)
  : light_environment(light_environment) {
}

EnvironmentEBT::EnvironmentEBT(spline::Spline *light_environment,
			       Disturbance *disturbance)
  : Environment(light_environment),
    disturbance(disturbance),
    age(0.0) {
}

double EnvironmentEBT::patch_survival(double age_at_birth) {
  return disturbance->survival_probability(age_at_birth, age);
}

}
