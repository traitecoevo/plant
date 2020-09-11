// Built from  inst/include/plant/models/ff16_strategy.h on Fri Jul 24 10:23:19 2020 using the scaffolder, from the strategy:  FF16
// -*-c++-*-
#ifndef PLANT_PLANT_K93_STRATEGY_H_
#define PLANT_PLANT_K93_STRATEGY_H_

#include <memory>
#include <plant/control.h>
#include <plant/qag_internals.h> // quadrature::intervals_type
#include <plant/internals.h> // quadrature::intervals_type
#include <plant/strategy.h>
#include <plant/models/k93_environment.h>

namespace plant {

class K93_Strategy: public Strategy<K93_Environment> {
public:
  typedef std::shared_ptr<K93_Strategy> ptr;
  K93_Strategy();

  // update this when the length of state_names changes
  static size_t state_size () { return 3; }
  // update this when the length of aux_names changes
  size_t aux_size () { return aux_names().size(); }

  static std::vector<std::string> state_names() {
    return  std::vector<std::string>({
      "height",
      "mortality",
      "fecundity"
      });
  }

  std::vector<std::string> aux_names() {
    return std::vector<std::string>({"competition_effect"});
  }

  void compute_rates(const K93_Environment& environment,
                     bool reuse_intervals,
                     Internals& vars);



  void refresh_indices();

  // These are stuck in plant.h
  double establishment_probability(const K93_Environment& environment);
  double net_mass_production_dt(const K93_Environment& environment,
                                double size, double cumulative_basal_area,
                                bool reuse_intervals=false);

  double Q(double z, double size) const;

  double compute_competition(double z, double size) const;

  void update_dependent_aux(const int index, Internals& vars);


  // K93 Methods  ----------------------------------------------
  double size_to_basal_area(double size) const;

  // Initial seedling size (dbh cm)
  double height_0;

  // * Growth
  // Growth intercept (yr-1)
  double b_0;
  // Growth asymptote (yr-1.(ln cm)-1)
  double b_1;
  // Growth suppression rate (m2.cm-2.yr-1)
  double b_2;

  // Growth rate of a plant per unit time:
  double size_dt(double size, double cumulative_basal_area) const;

  // * Reproduction
  // Recruitment rate (cm2.yr-1)
  double d_0;
  // Reduction from suppression (m2.cm-2.yr-1)
  double d_1;

  // Rate of offspring production
  double fecundity_dt(double size, double cumulative_basal_area) const;

  // * Mortality
  // Intercept (yr-1)
  double c_0;
  // Suppression rate (m2.cm-2.yr-1)
  double c_1;

  // Probability of survival during dispersal
  // required by scm.h
  double S_D = 1.0;

  // Smoothing parameter
  double eta = 12;

  // Rate of mortality over time
  double mortality_dt(double cumulative_basal_area,
                      double cumulative_mortality) const;

  // Set constants within K93_Strategy
  void prepare_strategy();

  std::string name;
};

K93_Strategy::ptr make_strategy_ptr(K93_Strategy s);

}

#endif
