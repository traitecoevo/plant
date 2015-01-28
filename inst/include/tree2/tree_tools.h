// -*-c++-*-
#ifndef TREE2_TREE_TOOLS_H_
#define TREE2_TREE_TOOLS_H_

namespace tree2 {
namespace tools {
Environment fixed_environment(double canopy_openness,
			      double height_max=150.0);
double lcp_whole_plant(Plant p);

}
}

#endif
