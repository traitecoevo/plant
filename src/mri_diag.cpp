// Diagnostic counter for the multirate (method="mri") fast sub-cycle.
//
// Each Patch::fast_rates call evaluates the soil (fast-block) coupling once,
// which costs O(M) cohort physiology solves -- the fast-block cost driver. This
// counter accumulates those calls across a run so we can read the per-macro-step
// fast-evaluation count (Lever-1 headroom: how many micro evaluations the stiff
// drainage currently forces) and compare schemes. Not part of the model; used
// only by the multirate cost measurements.
#include <Rcpp.h>

namespace plant {
long mri_fast_rate_calls = 0;
// Full patch RHS evaluations (Patch::compute_rates), i.e. one evaluation of the
// whole coupled derivative -- the cost/step proxy shared by every global stepper
// (rkck, rodas). Lets us compare accepted-step economics across methods and
// diagnose whether the coupled system is accuracy- or stability-limited.
long patch_rhs_calls = 0;
}

// [[Rcpp::export]]
double mri_fast_rate_calls_get() {
  return static_cast<double>(plant::mri_fast_rate_calls);
}

// [[Rcpp::export]]
void mri_fast_rate_calls_reset() {
  plant::mri_fast_rate_calls = 0;
}

// [[Rcpp::export]]
double patch_rhs_calls_get() {
  return static_cast<double>(plant::patch_rhs_calls);
}

// [[Rcpp::export]]
void patch_rhs_calls_reset() {
  plant::patch_rhs_calls = 0;
}
