
class Strategy {
public:
  typedef std::shared_ptr<Strategy> ptr;
  Strategy();

  // update this when the length of state_names changes
  static size_t state_size ();
  // update this when the length of aux_names changes
  size_t aux_size ();

  static std::vector<std::string> state_names();

  std::vector<std::string> aux_names();

  std::map<std::string, int> state_index; 
  std::map<std::string, int> aux_index; 
  
  bool collect_all_auxillary;

  void refresh_indices();

  double competition_effect(double size) const;
  double competition_effect_state(Internals& vars);

  void compute_rates(const Environment& environment, bool reuse_intervals,
                     Internals& vars);

  void update_dependent_aux(const int index, Internals& vars);

  double net_mass_production_dt(const Environment& environment,
                                double size, double competition_effect_,
                                bool reuse_intervals=false);

  double germination_probability(const Environment& environment);

  double compute_competition(double z, double size, double competition_effect) const;

  double initial_size(void) const;

  double size_0;

  // Every Strategy needs a set of Control objects
  Control control;

  std::string name;
};

Base_Strategy::ptr make_strategy_ptr(Base_Strategy s);

}


class Base_Strategy {
public:
  typedef std::shared_ptr<Base_Strategy> ptr;
  FF16_Strategy();

  // update this when the length of state_names changes
  static size_t state_size ();
  // update this when the length of aux_names changes
  size_t aux_size ();

  static std::vector<std::string> state_names();

  std::vector<std::string> aux_names();

  // TODO : expose this so can access state_names directly
  std::map<std::string, int> state_index; 
  std::map<std::string, int> aux_index; 
  
  bool collect_all_auxillary;

  // This would get renamed to something more generic - get_competiton_effect ?
  double area_leaf(double height) const;

  void refresh_indices();

  // This would get renamed to something more generic -- compute_rates
  void compute_vars_phys(const Environment& environment, bool reuse_intervals,
                Internals& vars);

  void update_dependent_aux(const int index, Internals& vars);


  // [eqn 20] Survival of seedlings during germination
  double germination_probability(const Environment& environment);

  // * Competitive environment
  // This would get renamed to something more generic -- compute_competition
  double area_leaf_above(double z, double height, double area_leaf_) const;

  // This would get renamed to something more generic -- initial size?
 
  // The aim is to find a plant height that gives the correct seed mass.
  double height_seed(void) const;

  // Every Strategy needs a set of Control objects
  Control control;

  std::string name;
};

Base_Strategy::ptr make_strategy_ptr(Base_Strategy s);

}
