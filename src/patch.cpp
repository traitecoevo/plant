#include "util.h"

#include "patch.h"

namespace model {

Patch::Patch(Parameters p)
  : standalone(true),
    parameters(new Parameters(p)) {
  set_strategies();
}

Patch::Patch(Parameters *p)
  : standalone(false),
    parameters(p) {
  set_strategies();
}

Patch::Patch(const Patch &other)
  : standalone(other.standalone) {
  Rprintf("Copy constructor\n");
  if ( standalone )
    parameters = new Parameters(*other.parameters);
  else
    parameters = other.parameters;
  set_strategies();
}

Patch& Patch::operator=(const Patch &rhs) {
  Rprintf("Assigmnent operator\n");
  // TODO: Violates DRY - must be some way of doing both.
  standalone = rhs.standalone;
  if ( standalone )
    parameters = new Parameters(*rhs.parameters);
  else
    parameters = rhs.parameters;
  set_strategies();

  return *this;
}

Patch::~Patch() {
  if ( standalone )
    delete parameters;
}

size_t Patch::size() const {
  return species.size();
}

Rcpp::List Patch::get_plants(int idx) const {
  util::check_bounds(idx, size());
  return species[idx].get_plants();
}

void Patch::add_seed(int idx) {
  species[idx].add_seed();
}

void Patch::r_add_seed(int idx) {
  util::check_bounds(idx, size());
  add_seed(idx);
}

void Patch::set_strategies() {
  species.clear();

  // This is really ugly, and I don't know that it is correct.
  for ( std::vector<Strategy>::iterator 
	  it = parameters->strategies.begin();
	it != parameters->strategies.end(); it++ ) {
    Species s(&(*it)); // ugly (iterator -> object -> pointer)
    species.push_back(s);
  }
    
}

}
