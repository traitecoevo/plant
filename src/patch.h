// -*-c++-*-
#ifndef TREE_PATCH_
#define TREE_PATCH_

#include <vector>

#include "parameters.h"
#include "species.h"

namespace model {

class Patch {
public:
  Patch(Parameters p);
  Patch(Parameters *p);

  Patch(const Patch &other);
  Patch& operator=(const Patch &rhs);
  ~Patch();

private:
  void set_strategies();

  bool standalone;
  Parameters *parameters;

  std::vector< Species > species;
};

}

#endif
