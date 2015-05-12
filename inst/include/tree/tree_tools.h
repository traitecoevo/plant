// -*-c++-*-
#ifndef TREE_TREE_TOOLS_H_
#define TREE_TREE_TOOLS_H_

namespace tree {
namespace tools {
Environment fixed_environment(double canopy_openness,
			      double height_max=150.0);
double lcp_whole_plant(FFW16_PlantPlus p);

}
}

#endif
