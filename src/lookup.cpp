#include "lookup.h"

namespace utils {

Lookup::Lookup() { 
  build_lookup();
}

Lookup::Lookup(const Lookup &other) {
  build_lookup();
  // not sure here -- this is slow but works.
  set_parameters(other.get_parameters());
  // might work better
  copy_parameters(other);
}

Rcpp::List Lookup::get_parameters() const {
  Rcpp::List ret;
  for ( lookup_type::const_iterator it = lookup_table.begin();
	it != lookup_table.end(); it++ )
    ret[it->first] = *(it->second);
  return ret;
}

void Lookup::set_parameters(Rcpp::List x) {
  if ( x.size() == 0 )
    return;

  std::vector<std::string> names = x.names();
  for ( int i = 0; i < x.size(); i++ ) {
    double *tmp = lookup(names[i]);
    *tmp = Rcpp::as<double>(x[i]);
  }
}

void Lookup::copy_parameters(const Lookup &source) {
  // Might also be possible to do this:
  // lookup_type::const_iterator it_from = source.lookup_table.begin();
  // lookup_type::iterator it_to = lookup_table.begin();
  // while ( it_to != lookup_table.end() ) {
  //   *(it_to->second) = *(it_from->second);
  //   it_to++;
  //   it_from++;
  // }
  for ( lookup_type::iterator it = lookup_table.begin();
	it != lookup_table.end(); it++ )
    *(it->second) = *(source.lookup(it->first));
}

void Lookup::build_lookup() {
  do_build_lookup();
}

double* Lookup::lookup(std::string key) const {
  lookup_type::const_iterator it = lookup_table.find(key);
  if ( it == lookup_table.end() )
    Rf_error("Key %s not found", key.c_str());
  return it->second;
}

}
