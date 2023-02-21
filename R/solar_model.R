# Returns incoming PAR on surface normal to incoming light, for given solar
# angle (in radians) S_atmos = radiation outside atmosphere W / m2/s tau:
# transmisivity of shortest atmospheric path PAR_frac Fraction of light that is
# Photosynthetically Active return value in mol PAR / m2/ day
PAR_given_solar_angle <- function(solar_angle, S_atmos = 1360, tau = 0.8, PAR_frac = 0.5) {
  
  # Adjust intensity for relative path-length through atmosphere

  # Would adjust for altitude here, but this function disabled
  altitude_adjustment <- 1
  
  M <- 1/cos(pi/2 - solar_angle) * altitude_adjustment
  
  mol_per_watt <- 4.6/1e+06
  sec_per_day <- 3600 * 24
  # Watts /m2/s
  ret <- (S_atmos * PAR_frac) * tau^M * mol_per_watt * sec_per_day
  ret[solar_angle <= 0] <- 0
  
  ret
}

# Returns solar angle = angle between the horizontal and the line connecting to
# the sun for given latitude and time of day, in radians
solar_angle <- function(decimal_day_time, latitude) {

  day <- floor(decimal_day_time)
  time <- (decimal_day_time - day) * 24
  
  radians <- function(x) {
    x/180 * pi
  }
  
  lat_radians <- radians(latitude)
  
  # solar declination (radians) - angle b/w the earth-sun line and the equatorial
  # plane
  delta <- 0.39785 * sin(4.869 + day/365 * 2 * pi + 0.03345 * sin(6.224 + day/365 * 
    2 * pi))
  
  # Hour Angle (degrees) - angular distance that earth has rotated in a day
  HrAng <- radians(180 - time * 15)
  
  # Sin of solar angle & Solar Angle - angle between the horizontal and the line
  # connecting to the sun
  SinA <- sin(lat_radians) * sin(delta) + cos(lat_radians) * cos(delta) * cos(HrAng)
  asin(SinA)
}

#' Title Calculate the hour of sun rise
#'
#' @param decimal_day_time Day of the year in decimal time
#' @param latitude Latitude in degrees
#'
#' @return Hour of sunrise (hr)
#' @export
#'
#' @examples sun_rise(90, -27)
#' 
sun_rise <- function(decimal_day_time, latitude){

  radians <- function(x) {
    x/180 * pi
  }
  lat_radians <- radians(latitude)
  day <- floor(decimal_day_time)
  delta <- 0.39785 * sin(4.869 + day/365 * 2 * pi + 0.03345 * sin(6.224 + day/365 * 
                                                                    2 * pi))
  
  Hrang <- acos(-1*(sin(lat_radians)*sin(delta)) / (cos(lat_radians)*cos(delta)))
  sunrise <- (Hrang / pi * 180 - 180) / -15
  sunrise
}
