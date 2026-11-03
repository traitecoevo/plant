// -*-c++-*-
#ifndef PLANT_PLANT_TF24_DIURNAL_H_
#define PLANT_PLANT_TF24_DIURNAL_H_

#include <plant/util.h>
#include <cmath>
#include <string>

namespace plant {
namespace diurnal {

// The shape of a day, as TF24's daily integral reads it (#618, experiment).
//
// Time is patch age in years with t = 0 at 1 January, and the day of year is
// 365 * frac(t). It is kept CONTINUOUS rather than floored to a whole day, so
// every quantity here is smooth in t and an adaptive step never meets a jump at
// midnight. Hours are local solar time, so noon is 12 by definition and the day
// is symmetric about it.

constexpr double pi = 3.14159265358979323846;

// Solar declination (radians).
//
// ⚠️ THE EXPRESSION IS sin(declination), NOT the declination: its amplitude
// 0.39785 is sin(23.44 deg), the obliquity. plant's R solar_angle() uses the
// same expression as the angle itself, which understates the solstice
// declination by about 0.6 degrees; this takes the arcsine.
inline double declination(double time) {
  const double x = 2.0 * pi * (time - std::floor(time));
  return std::asin(0.39785 * std::sin(4.869 + x + 0.03345 * std::sin(6.224 + x)));
}

// The sun on one day at one latitude, reduced to what the integral reads:
// sin(elevation) at hour h is a + b cos(hour angle), and the photoperiod is
// where that is positive.
struct SolarDay {
  double a = 0.0;          // sin(latitude) sin(declination)
  double b = 0.0;          // cos(latitude) cos(declination)
  double daylength = 0.0;  // hours from sunrise to sunset

  SolarDay(double time, double latitude) {
    if (!util::is_finite(latitude) || std::abs(latitude) > 90.0) {
      util::stop("diurnal: latitude must be a finite number of degrees in "
                 "[-90, 90]; got " + util::to_string(latitude));
    }
    const double phi = latitude * pi / 180.0;
    const double delta = declination(time);
    a = std::sin(phi) * std::sin(delta);
    b = std::cos(phi) * std::cos(delta);
    // cos(hour angle at sunrise) = -tan(latitude) tan(declination). Outside
    // (-1, 1) the sun never rises or never sets, and there is no photoperiod to
    // integrate over.
    const double cos_sunrise = -a / b;
    if (!(cos_sunrise > -1.0 && cos_sunrise < 1.0)) {
      util::stop("diurnal: no sunrise or no sunset at latitude " +
                 util::to_string(latitude) + " on day " +
                 util::to_string(365.0 * (time - std::floor(time))) +
                 " (polar day or night)");
    }
    daylength = 24.0 / pi * std::acos(cos_sunrise);
  }

  double sin_elevation(double hour) const {
    return a + b * std::cos(pi * (hour - 12.0) / 12.0);
  }
};

// Relative diurnal course of the vapour pressure deficit: one cosine harmonic
// peaking at `peak_hour`. Divided by its value at noon it is the multiplier on
// the noon VPD, so noon reads the driver exactly.
inline double vpd_shape(double hour, double amplitude, double peak_hour) {
  return 1.0 + amplitude * std::cos(2.0 * pi * (hour - peak_hour) / 24.0);
}

}
}

#endif
