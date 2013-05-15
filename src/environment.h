// -*-c++-*-
#ifndef TREE_ENVIRONMENT_H_
#define TREE_ENVIRONMENT_H_

#include "spline.h"
#include "disturbance.h"

namespace model {

class Environment {
public:
  Environment(spline::Spline *light_environment);
  virtual ~Environment() {}
  spline::Spline *light_environment;
};

class EnvironmentEBT : public Environment {
  EnvironmentEBT(spline::Spline *light_environment,
		 Disturbance *disturbance);
  double patch_survival(double age_at_birth);
private:
  Disturbance *disturbance;
  double age;
};

}

#endif
