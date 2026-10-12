// -*-c++-*-
#ifndef PLANT_PLANT_WINDOW_REGIME_H_
#define PLANT_PLANT_WINDOW_REGIME_H_

#include <vector>
#include <algorithm>
#include <plant/disturbance_regime.h>
#include <plant/util.h>

namespace plant {

// The weighting a continuous patch gives an arrival (#519). Where the Weibull
// regime weights a cohort by the density of patches at its age, this weights it
// by WHEN it arrived in one long run: a tabulated weight over an arrival window
// [start, end], zero outside it. The patch normalises it over the cohorts the
// run introduced (is_arrival_weighting), so the lifetime-fitness integral over
// birth dates reads as the weighted mean of a cohort's lifetime output across
// the window's arrivals, exactly, however the introduction grid sits against
// the window's edges; and a resident's equilibrium is the same fixed point it
// is under the Weibull. The normalisation here, over the window itself, only
// keeps density() on the scale of 1 / (end - start) for a reader.
//
// No knots means a uniform window, density 1 / (end - start). Knots give a
// piecewise-linear weight, read at the cohort's own birth date and held at the
// end knots' values beyond them: a germination opportunity built from the same
// rainfall series that drives the run, say. The regime never reads the
// environment; what it weights by is the caller's to tabulate.
//
// Survival is one: nothing resets the patch, so a cohort's own mortality is the
// only hazard it meets.
class Window_Regime : public Disturbance_Regime {
public:
  Window_Regime(double start, double end,
                std::vector<double> x, std::vector<double> w)
    : start_(start), end_(end), x_(std::move(x)), w_(std::move(w)) {
    if (!(end_ > start_)) {
      util::stop("window_end must exceed window_start");
    }
    if (x_.size() != w_.size()) {
      util::stop("window_weight_x and window_weight_y must have the same length");
    }
    if (!util::is_sorted(x_.begin(), x_.end())) {
      util::stop("window_weight_x must be increasing");
    }
    for (const double v : w_) {
      if (!(v >= 0.0)) {
        util::stop("window weights must be non-negative");
      }
    }
    // The integral of the weight over the window, by the trapezium over the
    // knots inside it and its two ends, which is exact for a piecewise-linear
    // weight.
    std::vector<double> knots{start_};
    for (const double t : x_) {
      if (t > start_ && t < end_) {
        knots.push_back(t);
      }
    }
    knots.push_back(end_);
    double total = 0.0;
    for (size_t i = 1; i < knots.size(); ++i) {
      total += (knots[i] - knots[i - 1]) *
               (weight(knots[i - 1]) + weight(knots[i])) / 2.0;
    }
    if (!(total > 0.0)) {
      util::stop("the window weights integrate to zero over the window");
    }
    norm_ = 1.0 / total;
  }

  virtual double density(double time) const {
    return (time >= start_ && time <= end_) ? weight(time) * norm_ : 0.0;
  }
  virtual double pr_survival(double) const { return 1.0; }
  virtual bool is_arrival_weighting() const { return true; }

private:
  // The tabulated weight at t, before normalisation.
  double weight(double t) const {
    if (x_.empty()) {
      return 1.0;
    }
    if (t <= x_.front()) {
      return w_.front();
    }
    if (t >= x_.back()) {
      return w_.back();
    }
    // The first knot strictly above t; the one below it is at or below t, and
    // strictly below the one above, so the width is never zero.
    const size_t i = std::upper_bound(x_.begin(), x_.end(), t) - x_.begin();
    const double f = (t - x_[i - 1]) / (x_[i] - x_[i - 1]);
    return w_[i - 1] + f * (w_[i] - w_[i - 1]);
  }

  double start_;
  double end_;
  double norm_;
  std::vector<double> x_;
  std::vector<double> w_;
};

}

#endif
