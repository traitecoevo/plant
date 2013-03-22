// -*-c++-*-
#ifndef TREE_POPULATION_
#define TREE_POPULATION_

#include "parameters.h"

namespace model {

class Population {
public:
  Population(Parameters p);
  Population(Parameters *p);

  Population(const Population &other);
  Population& operator=(const Population &rhs);
  ~Population();

private:
  bool standalone;
  Parameters *parameters;
};

}

#endif
