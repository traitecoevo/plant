#include "lookup.h"

namespace utils {

Lookup::Lookup() : addr(NULL) { }

Rcpp::List Lookup::get_parameters() {
  build_lookup();
  Rcpp::List ret;
  for ( lookup_type::const_iterator it = lookup_table.begin();
	it != lookup_table.end(); it++ )
    ret[it->first] = *(it->second);
  return ret;
}

void Lookup::set_parameters(Rcpp::List x) {
  if ( x.size() == 0 )
    return;
  build_lookup();
  std::vector<std::string> names = x.names();
  for ( int i = 0; i < x.size(); i++ ) {
    double *tmp = lookup(names[i]);
    *tmp = Rcpp::as<double>(x[i]);
  }
}

void Lookup::build_lookup() {
  // The first condition detects copies, the second detects that the
  // list was never built.
  if ( addr != this ) {
    do_build_lookup();
    addr = this;
  }
}

double* Lookup::lookup(std::string key) const {
  lookup_type::const_iterator it = lookup_table.find(key);
  if ( it == lookup_table.end() )
    Rf_error("Key %s not found", key.c_str());
  return it->second;
}

}
