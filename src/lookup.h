// -*-c++-*-
#ifndef TREE_LOOKUP_
#define TREE_LOOKUP_

#include <string>
#include <map>
#include <Rcpp.h>

// This is a weird class, partly because the way it is used is weird.
// 
// We have a class with an absolute ton of "double" members
// (`Strategy`); these represent model parameters and we want to be
// able to access them easily from within the `Plant` class (`Plant`
// has a member `Strategy *strategy` so we can do `strategy->n_area`
// to get concentration of leaf nitrogen by area).
// 
// However, we need to be able to get and set these parameters from R,
// but I don't want to write a thousand small interfaces -- even with
// Rcpp modules that's going to be a pain.  I want to be able to do
// get_parameters() and have a list of parameters (in this case, an R
// list) and similarly do set_parameters(x) for some list x and set a
// subset of parameters. 
// 
// That's what this class does.  Inherit from it privately and define
// the method `do_build_lookup()`.  This ends up with lines that look
// like this:
//   lookup_table["foo"] = &foo;
// where "foo" is the name of the element in the R list and &foo is
// the address of the member (which *must* be a double).

// TODO: It might be easiest to do the build_lookup at get/set
// parameter time (which is never time sensitive) than at copy (which
// might conceivably be).

namespace utils {

class Lookup {
  Lookup();
  Lookup(const Lookup &other);
  // Assignment operator?
  Rcpp::List get_parameters() const;
  void set_parameters(Rcpp::List x);

private:
  typedef std::map<std::string, double*> lookup_type;
  lookup_type lookup_table;
  void build_lookup(); // used on construction
  virtual void do_build_lookup() = 0;
  double* lookup(std::string key) const;
  void copy_parameters(const Lookup &source);
};

}

#endif
