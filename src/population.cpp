#include "population.h"

namespace model {

Population::Population(Parameters p)
  : standalone(true),
    parameters(new Parameters(p)) {
}

Population::Population(Parameters *p)
  : standalone(false),
    parameters(p) { 
}

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
  // TODO: Violates DRY - must be some way of doing both.
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
