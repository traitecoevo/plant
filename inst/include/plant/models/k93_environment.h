// -*-c++-*-
#ifndef PLANT_PLANT_K93_ENVIRONMENT_H_
#define PLANT_PLANT_K93_ENVIRONMENT_H_

#include <plant/models/ff16_environment.h>

using namespace Rcpp;

namespace plant {

class K93_Environment : public FF16_Environment {
public:

  // Reuse constructors from FF16_Environment
  using FF16_Environment::FF16_Environment;

};

}

#endif
