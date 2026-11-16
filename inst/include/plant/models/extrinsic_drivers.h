#ifndef PLANT_EXTRINSIC_DRIVERS_H
#define PLANT_EXTRINSIC_DRIVERS_H

#include <odelia/drivers.hpp>

namespace plant {

// The drivers a strategy reads: a named constant, or a series supplied from
// outside the model and read between its control points.
//
// One type, and it lives beside the interpolant it reads through. There was a
// second copy of it here, class for class and statement for statement, which could
// only ever agree with this one by being edited twice.
using ExtrinsicDrivers = odelia::drivers::Drivers;

}

#endif //PLANT_EXTRINSIC_DRIVERS_H
