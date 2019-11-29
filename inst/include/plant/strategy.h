// -*-c++-*-
#ifndef PLANT_PLANT_STRATEGY_H_
#define PLANT_PLANT_STRATEGY_H_

#include <memory>
#include <plant/control.h>
#include <plant/qag_internals.h> // quadrature::intervals_type
#include <plant/internals.h> // quadrature::intervals_type

namespace plant {

// Environment needs Parameters to initialise and that needs Strategy,
// so there's a really awkward circular reference here.  This forward
// declaration breaks it, but there might be a better solution.
class Environment;

class Strategy {
public:

  /* Common variables */

  typedef std::shared_ptr<Strategy> ptr;

  std::map<std::string, int> state_index; 
  std::map<std::string, int> aux_index; 

  bool collect_all_auxillary;

  double size_0;

  // Control objects govern things to do with how numerical calculations 
  // are performed, rather than the biological control that this class has.
  Control control;

  std::string name;


  /* Interface methods that all strategies must define. */

  void compute_rates(const Environment& environment, bool reuse_intervals,
                Internals& vars);

  double compute_competition(double z, double size, double competition_effect) const;

  virtual double competition_effect(double size) const = 0;

  double net_mass_production_dt(const Environment& environment,
                                double size, double competition_effect_,
                                bool reuse_intervals=false);

  double germination_probability(const Environment& environment);

  double fecundity_dt(double net_mass_production_dt,
                      double fraction_allocation_reproduction) const;

  double mortality_dt(double productivity_area, double cumulative_mortality) const;


  /* State names. Generally overridden. Must be a superset containing these in the same order */

  static std::vector<std::string> state_names() {
    return  std::vector<std::string>({
      "height",
      "mortality",
      "fecundity",
    });
  }

  std::vector<std::string> aux_names() {
    return std::vector<std::string>({
      "competition_effect",
      "net_mass_production_dt"
    });
  }


  /* General utility methods */

  // update this when the length of state_names changes
  static size_t state_size ();
  // update this when the length of aux_names changes
  size_t aux_size () { return aux_names().size(); }


  void refresh_indices () {
      // Create and fill the name to state index maps
    state_index = std::map<std::string,int>();
    aux_index   = std::map<std::string,int>();
    std::vector<std::string> aux_names_vec = aux_names();
    std::vector<std::string> state_names_vec = state_names();
    for (int i = 0; i < state_names_vec.size(); i++) {
      state_index[state_names_vec[i]] = i;
    }
    for (int i = 0; i < aux_names_vec.size(); i++) {
      aux_index[aux_names_vec[i]] = i;
    }
  }

  // for updating auxillary state
  void update_dependent_aux(const int index, Internals& vars) {
    if (index == HEIGHT_INDEX) {
      double height = vars.state(HEIGHT_INDEX);
      vars.set_aux(aux_index.at("competition_effect"), competition_effect(height));
    }
  }

};

}

#endif
