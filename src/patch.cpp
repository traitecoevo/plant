#include "patch.h"

namespace model {

Patch::Patch(Parameters p)
  : standalone(true),
    parameters(new Parameters(p)) {
}

Patch::Patch(Parameters *p)
  : standalone(false),
    parameters(p) { 
}

Patch::Patch(const Patch &other)
  : standalone(other.standalone) {
  Rprintf("Copy constructor\n");
  if ( standalone )
    parameters = new Parameters(*other.parameters);
  else
    parameters = other.parameters;
}

Patch& Patch::operator=(const Patch &rhs) {
  Rprintf("Assigmnent operator\n");
  // TODO: Violates DRY - must be some way of doing both.
  standalone = rhs.standalone;
  if ( standalone )
    parameters = new Parameters(*rhs.parameters);
  else
    parameters = rhs.parameters;

  return *this;
}

Patch::~Patch() {
  if ( standalone )
    delete parameters;
}

}
