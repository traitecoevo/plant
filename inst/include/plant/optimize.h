// -*-c++-*-
#ifndef PLANT_PLANT_OPTIMIZE_H_
#define PLANT_PLANT_OPTIMIZE_H_

// Robust 1-D function minimiser: Brent's method (golden-section search with
// parabolic interpolation). This is the algorithm behind R's optimize() /
// Forsythe-Malcolm-Moler `fmin`. It converges super-linearly on a smooth
// objective near its optimum but falls back to a golden-section step whenever
// the parabolic fit is unreliable, so it keeps the bracketing robustness that
// the leaf hydraulic solvers depend on (cf. the TOMS748 caveat in uniroot.h:
// here we only ever evaluate STRICTLY INTERIOR points, never the clamped
// bracket endpoints).

#include <cfloat>  // DBL_EPSILON
#include <cmath>

namespace plant {
namespace util {

// Minimise f over [ax, bx] (requires ax <= bx). Returns the argmin; on return
// fmin (if non-null) holds f at the argmin. `tol` is the absolute tolerance on
// the location of the minimum; values below ~sqrt(DBL_EPSILON)*|x| are not
// useful. A direct transcription of R's Brent_fmin (src/appl/fmin.c).
template <typename Function>
double brent_fmin(Function f, double ax, double bx, double tol,
                  double* fmin = nullptr) {
  // c is the squared inverse of the golden ratio (~0.3819660).
  const double c = (3.0 - std::sqrt(5.0)) * 0.5;
  const double eps = std::sqrt(DBL_EPSILON);

  double a = ax, b = bx;
  double v = a + c * (b - a);
  double w = v, x = v;
  double d = 0.0, e = 0.0;
  double fx = f(x);
  double fv = fx, fw = fx;
  const double tol3 = tol / 3.0;

  for (;;) {
    const double xm = (a + b) * 0.5;
    const double tol1 = eps * std::abs(x) + tol3;
    const double tol2 = tol1 * 2.0;

    // Convergence check.
    if (std::abs(x - xm) <= tol2 - (b - a) * 0.5)
      break;

    double p = 0.0, q = 0.0, r = 0.0;
    bool use_golden = true;

    if (std::abs(e) > tol1) {
      // Fit a parabola through (x,fx), (v,fv), (w,fw).
      r = (x - w) * (fx - fv);
      q = (x - v) * (fx - fw);
      p = (x - v) * q - (x - w) * r;
      q = (q - r) * 2.0;
      if (q > 0.0)
        p = -p;
      else
        q = -q;
      r = e;
      e = d;
      // Accept the parabolic step only if it is well inside (a,b) and shrinks.
      if (std::abs(p) < std::abs(0.5 * q * r) &&
          p > q * (a - x) && p < q * (b - x)) {
        d = p / q;
        const double u = x + d;
        // Keep the new point away from the bracket endpoints.
        if (u - a < tol2 || b - u < tol2)
          d = (x < xm) ? tol1 : -tol1;
        use_golden = false;
      }
    }

    if (use_golden) {
      e = (x < xm) ? (b - x) : (a - x);
      d = c * e;
    }

    // Evaluate f at a point at least tol1 away from x.
    double u;
    if (std::abs(d) >= tol1)
      u = x + d;
    else
      u = (d > 0.0) ? (x + tol1) : (x - tol1);
    const double fu = f(u);

    // Update the bracket and the three best points.
    if (fu <= fx) {
      if (u < x) b = x; else a = x;
      v = w; fv = fw;
      w = x; fw = fx;
      x = u; fx = fu;
    } else {
      if (u < x) a = u; else b = u;
      if (fu <= fw || w == x) {
        v = w; fv = fw;
        w = u; fw = fu;
      } else if (fu <= fv || v == x || v == w) {
        v = u; fv = fu;
      }
    }
  }

  if (fmin != nullptr)
    *fmin = fx;
  return x;
}

// Golden-section search for the MAXIMUM of a unimodal f over [ax, bx]. Returns
// the argmax (midpoint of the final bracket); terminates when the bracket width
// falls to `tol`. Reuses one interior golden point per iteration, so it costs a
// single new f() evaluation per step after the initial two.
//
// Why this exists alongside brent_fmin: brent_fmin converges faster but its
// parabolic step makes the argmax a *non-smooth* function of the inputs. The
// production collar solver (Leaf::find_root_collar_psi) feeds its argmax into the
// demographic growth-rate gradient, which needs that argmax to vary smoothly
// with plant state -- so it uses this fixed-iteration golden-section search.
// Where the argmax does not feed a gradient (the single-layer leaf optimisers),
// prefer brent_fmin, which was measured ~2.3-2.6x faster there.
template <typename Function>
double golden_section_max(Function f, double ax, double bx, double tol) {
  const double gr = (std::sqrt(5.0) + 1.0) / 2.0;  // ~1.6180339...
  double a = ax, b = bx;
  double c = b - (b - a) / gr;
  double d = a + (b - a) / gr;
  double fc = f(c);
  double fd = f(d);
  while (std::abs(b - a) > tol) {
    if (fc > fd) {
      b  = d;
      d  = c;
      fd = fc;                 // reuse
      c  = b - (b - a) / gr;
      fc = f(c);               // 1 new eval
    } else {
      a  = c;
      c  = d;
      fc = fd;                 // reuse
      d  = a + (b - a) / gr;
      fd = f(d);               // 1 new eval
    }
  }
  return (a + b) / 2.0;
}

}
}

#endif
