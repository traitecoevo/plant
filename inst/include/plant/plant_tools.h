// -*-c++-*-
#ifndef PLANT_PLANT_TOOLS_H_
#define PLANT_PLANT_TOOLS_H_

namespace plant {
namespace tools {
Environment fixed_environment(double canopy_openness,
			      double height_max=150.0);
double lcp_whole_plant(FFW16_PlantPlus p);

}
}

#endif
