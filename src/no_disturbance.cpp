#include <plant/disturbances/no_disturbance.h>

namespace plant {

double No_Disturbance::density(double time) const {
  return 1.0;
}

double No_Disturbance::pr_survival(double time) const {
  return 1.0;
}

}
