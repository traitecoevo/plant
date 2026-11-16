// -*-c++-*-
#ifndef PLANT_PLANT_DISTURBANCE_REGIME_H_
#define PLANT_PLANT_DISTURBANCE_REGIME_H_

#include <vector>
#include <RcppCommon.h> // NA_REAL

namespace plant {

class Disturbance_Regime {
public:
  virtual double density(double time) const {return NA_REAL;}
  virtual double pr_survival(double time) const {return NA_REAL;}
  virtual double cdf(double time) const {return NA_REAL;}
  virtual double icdf(double p) const {return NA_REAL;}
  virtual double r_mean_interval() const {return NA_REAL;}
  // Whether density() is a weighting of THIS run's arrivals, to be normalised
  // over the cohorts the run actually introduced, rather than a density over
  // patch ages that integrates to one on its own. A weighting that steps to
  // zero at a window's edge would otherwise hand each edge cohort a whole
  // trapezium interval instead of half of one, and read the window mean a
  // factor 1 + spacing / width too high.
  virtual bool is_arrival_weighting() const {return false;}
  /* base class should define virtual destructor to avoid undefined behaviour 
  when static and dynamic types differ for an object (ie through pointers) and 
  an attempt to delete it is made */
  virtual ~Disturbance_Regime() = default;
  double r_density(double time) const {return density(time);}
};

}

#endif
