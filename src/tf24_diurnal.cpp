#include <plant.h>
#include <plant/models/tf24_diurnal.h>

// The sun as TF24's daily integral sees it (#618, experiment), for a script
// that builds noon drivers and must place them on the same day the model does.
// One row per time: the photoperiod (hours) and sin(elevation) at solar noon.
// [[Rcpp::export]]
Rcpp::DataFrame tf24_solar_day(std::vector<double> time, double latitude) {
  std::vector<double> daylength(time.size()), sin_elevation_noon(time.size());
  for (size_t i = 0; i < time.size(); ++i) {
    const plant::diurnal::SolarDay day(time[i], latitude);
    daylength[i] = day.daylength;
    sin_elevation_noon[i] = day.sin_elevation(12.0);
  }
  return Rcpp::DataFrame::create(Rcpp::Named("time") = time,
                                 Rcpp::Named("daylength") = daylength,
                                 Rcpp::Named("sin_elevation_noon") =
                                   sin_elevation_noon);
}

// sin(elevation) at each solar hour of one day.
// [[Rcpp::export]]
std::vector<double> tf24_sin_elevation(double time, std::vector<double> hour,
                                       double latitude) {
  const plant::diurnal::SolarDay day(time, latitude);
  std::vector<double> out(hour.size());
  for (size_t i = 0; i < hour.size(); ++i) {
    out[i] = day.sin_elevation(hour[i]);
  }
  return out;
}
