// -*-c++-*-
#ifndef TREE_PATCH_
#define TREE_PATCH_

#include "parameters.h"

namespace model {

class Patch {
public:
  Patch(Parameters p);
  Patch(Parameters *p);

  Patch(const Patch &other);
  Patch& operator=(const Patch &rhs);
  ~Patch();

private:
  bool standalone;
  Parameters *parameters;
};

}

#endif
