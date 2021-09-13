#include <plant/models/ff16_environment.h>

namespace plant {

FF16_Environment::FF16_Environment() {
  time = 0.0;

  canopy_light_tol = 1e-6;
  canopy_light_nbase = 17;
  canopy_light_max_depth = 16;
  canopy_rescale_usually = false;

  canopy = Canopy(canopy_light_tol, canopy_light_nbase, canopy_light_max_depth);

  soil_number_of_depths = 0;
  soil_initial_state = 0.0;
  soil_infiltration_rate = 1.0;

  vars = Internals(soil_number_of_depths);
  set_soil_water_state(soil_initial_state);
}

} // namespace plant
