#include "lookup.h"

namespace util {

Lookup::Lookup() : addr(NULL) { }

Rcpp::List Lookup::get_parameters() {
  build_lookup();
  Rcpp::List ret;
  for (lookup_type::const_iterator it = lookup_table.begin();
       it != lookup_table.end(); ++it)
    ret[it->first] = *(it->second);
  return ret;
}

void Lookup::set_parameters(Rcpp::List x) {
  if (x.size() == 0)
    return;
  build_lookup();
  std::vector<std::string> names = x.names();
  check_keys(names);

  for (int i = 0; i < x.size(); ++i) {
    double *tmp = lookup(names[static_cast<size_t>(i)]);
    *tmp = Rcpp::as<double>(x[i]);
  }
  set_parameters_post_hook();
}

bool Lookup::has_key(std::string key) const {
  lookup_type::const_iterator it = lookup_table.find(key);
  return it != lookup_table.end();
}

// Default hook is to do nothing.
void Lookup::set_parameters_post_hook() {
}

void Lookup::build_lookup() {
  if (addr != this) {
    do_build_lookup();
    addr = this;
  }
}

double* Lookup::lookup(std::string key) const {
  lookup_type::const_iterator it = lookup_table.find(key);
  if (it == lookup_table.end())
    ::Rf_error("Key %s not found", key.c_str());
  return it->second;
}

void Lookup::check_keys(std::vector<std::string> keys) const {
  std::vector<std::string>::const_iterator key = keys.begin();
  while (key != keys.end()) {
    if (!has_key(*key))
      ::Rf_error("Key %s not found", key->c_str());
    ++key;
  }
}

}
