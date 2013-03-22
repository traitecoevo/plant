#include "population.h"

namespace model {

Population::Population(Parameters p)
  : standalone(true),
    parameters(new Parameters(p)) {
  // TODO: What about the constants here?
}

Population::Population(Parameters *p)
  : standalone(false),
    parameters(p) { }

Population::Population(const Population &other)
  : standalone(other.standalone) {
  Rprintf("Copy constructor\n");
  if ( standalone )
    parameters = new Parameters(*other.parameters);
  else
    parameters = other.parameters;
}

Population& Population::operator=(const Population &rhs) {
  Rprintf("Assigmnent operator\n");
  // TODO: Violates DRY - must be some way of doing both.  This will
  // get more important once the state attributes have been added.
  standalone = rhs.standalone;
  if ( standalone )
    parameters = new Parameters(*rhs.parameters);
  else
    parameters = rhs.parameters;

  return *this;
}

Population::~Population() {
  if ( standalone )
    delete parameters;
}

}
