#include <plant/models/ff16_environment.h>

namespace plant {

  FF16_Environment::FF16_Environment(double canopy_light_tol = 1e-4,
                                     int canopy_light_nbase = 17,
                                     int canopy_light_max_depth = 16,
                                     bool canopy_rescale_usually = true,
                                     int soil_number_of_depths = 0,
                                     std::vector<double> soil_initial_state = std::vector<double>(1, 0.0),
                                     double infil_rate = 0.0) {
    time = 0.0;
    canopy = Canopy(canopy_light_tol, canopy_light_nbase, canopy_light_max_depth);
    vars = Internals(soil_number_of_depths);
    soil_infiltration_rate = infil_rate;
    set_soil_water_state(soil_initial_state);
  }

} // namespace plant
