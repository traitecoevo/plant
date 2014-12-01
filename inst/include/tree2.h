// -*-c++-*-
#ifndef _TREE2_H_
#define _TREE2_H_

#include <tree2/ode_bits.h>

// Include this early on.  It can be either after classes have been
// declared (but before Rcpp has been loaded) or first.  This file will
// attempt to provide declarations for the classes and namespaces that
// you use, but this might be fragile.
#include <tree2/RcppR6_pre.hpp>

// Extra shit:
namespace Rcpp {
SEXP wrap(const ode::state_saver<std::vector<double> >&);
}

// Anything after this point is OK to include Rcpp.h.  This is
// probably where the meat of the included material goes if your
// classes directly use Rcpp types.  Otherwise you can just declare
// them earlier up.

#include <Rcpp.h>
#include <tree2/lorenz.h>

// This line can safely be the last line in the file, but may go any
// point after RcppR6_pre.hpp is included.
#include <tree2/RcppR6_post.hpp>
#include <tree2/util_rcpp.h>

namespace Rcpp {
inline SEXP wrap(const ode::state_saver<std::vector<double> >& x) {
  Rcpp::NumericMatrix m(util::to_rcpp_matrix(x.y));
  m.attr("t") = x.t;
  return m;
}
}

#endif
