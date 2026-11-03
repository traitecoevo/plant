# TF24's daily integral (#618, experiment branch): the fixed twelve-hour day it
# always ran, one solve at solar noon scaled by c, or a diurnal integral.

tf24_day <- function(method, n = 48L, c_A = 1, c_E = 1, time = 0.2,
                     latitude = -33.6, amp = 0.3, height = 5, light = 1,
                     theta = 0.3) {
  s <- TF24_Strategy()
  s$control$daily_integration <- method
  s$control$daily_integration_points <- as.integer(n)
  s$control$one_point_c_assimilation <- c_A
  s$control$one_point_c_transpiration <- c_E
  env <- Environment("TF24")
  env$set_fixed_environment(light, height_max = 150)
  env$set_soil_water_state(rep(theta, env$get_soil_number_of_depths()))
  env$latitude <- latitude
  env$extrinsic_drivers_set_constant("atm_vpd_diurnal_amp", amp)
  env$time <- time
  ind <- TF24_Individual(s)
  ind$set_state("height", height)
  ind$compute_rates(env)
  list(P = ind$aux("net_mass_production_dt"), profit = ind$aux("profit"),
       transpiration = ind$aux("transpiration"),
       consumption = ind$internals$consumption_rates)
}

test_that("the sun is placed where the geometry says", {
  ## At the equator every day is twelve hours, and a latitude and its mirror
  ## share a day: acos(x) + acos(-x) = pi.
  expect_equal(tf24_solar_day(c(0, 0.3, 0.7), 0)$daylength, rep(12, 3),
               tolerance = 1e-14)
  t <- seq(0, 1, length.out = 13)
  expect_equal(tf24_solar_day(t, -33.6)$daylength +
                 tf24_solar_day(t, 33.6)$daylength, rep(24, 13),
               tolerance = 1e-14)
  ## Southern summer days are the long ones.
  expect_gt(tf24_solar_day(0, -33.6)$daylength,
            tf24_solar_day(0.5, -33.6)$daylength)
  ## The sun is on the horizon at sunrise and sunset.
  d <- tf24_solar_day(0.2, -33.6)$daylength
  expect_equal(tf24_sin_elevation(0.2, 12 + c(-1, 1) * d / 2, -33.6), c(0, 0),
               tolerance = 1e-12)
  expect_error(tf24_solar_day(0.5, 80), "polar")
})

test_that("one diurnal point is the one-point method with c = 1", {
  ## n = 1 places its single node at solar noon, so it is the noon solve billed
  ## over the photoperiod -- exactly.
  expect_identical(tf24_day("diurnal", n = 1L), tf24_day("one-point"))
})

test_that("the one-point constants scale carbon and water and nothing else", {
  base <- tf24_day("one-point")
  carbon <- tf24_day("one-point", c_A = 1.3)
  water <- tf24_day("one-point", c_E = 0.7)
  expect_identical(carbon$profit, 1.3 * base$profit)
  expect_identical(carbon$consumption, base$consumption)
  expect_identical(water$profit, base$profit)
  expect_equal(water$consumption, 0.7 * base$consumption, tolerance = 1e-14)
})

test_that("one-point and fixed-12h differ only by the photoperiod", {
  ## The same solve at the same noon drivers; one bills the day's photoperiod
  ## and the other twelve hours.
  one <- tf24_day("one-point")
  fixed <- tf24_day("fixed-12h")
  d <- tf24_solar_day(0.2, -33.6)$daylength
  expect_identical(one$profit, fixed$profit)
  expect_equal(one$consumption / fixed$consumption,
               rep(d / 12, length(one$consumption)), tolerance = 1e-14)
})

test_that("the diurnal integral has converged at 48 points", {
  p48 <- tf24_day("diurnal", n = 48L)$profit
  p96 <- tf24_day("diurnal", n = 96L)$profit
  expect_lt(abs(p48 - p96) / abs(p96), 1e-3)
  ## And a day is not its noon: with light falling away from noon, the daylight
  ## mean is below the noon rate.
  expect_lt(p96, tf24_day("one-point")$profit)
})

test_that("what the new methods cannot honour is refused", {
  expect_error(tf24_day("hourly"), "Unknown daily_integration")
  expect_error(tf24_day("diurnal", n = 0L), "at least 1")
  expect_error(tf24_day("one-point", c_A = 0), "must be positive")
  expect_error(tf24_day("diurnal", latitude = NaN), "latitude")
  expect_error(tf24_day("diurnal", amp = 1), "atm_vpd_diurnal_amp")
  ## The default never reads the site, so it runs without one.
  expect_silent(tf24_day("fixed-12h", latitude = NaN))
})

test_that("a light-dependent c reads the crown's noon light, clamped", {
  ## c(I) = c + s log(I'/500), I' the noon intercepted PPFD held to the
  ## calibrated range. At full openness I is k_I * PPFD_above.
  day <- function(...) {
    s <- TF24_Strategy()
    ctrl <- list(...)
    for (k in names(ctrl)) s$control[[k]] <- ctrl[[k]]
    env <- Environment("TF24")
    env$latitude <- -33.6
    env$time <- 0.2
    env$set_fixed_environment(1, height_max = 150)
    env$set_soil_water_state(rep(0.3, env$get_soil_number_of_depths()))
    ind <- TF24_Individual(s)
    ind$set_state("height", 5)
    ind$compute_rates(env)
    list(profit = ind$aux("profit"), k_I = s$pars$k_I,
         I = s$pars$k_I * 4.57 * env$f_PAR * env$get_SWdown())
  }
  base <- day(daily_integration = "one-point")
  lit <- day(daily_integration = "one-point", one_point_c_assimilation = 0.7,
             one_point_c_assimilation_slope = 0.15)
  expect_equal(lit$profit / base$profit, 0.7 + 0.15 * log(base$I / 500),
               tolerance = 1e-12)
  ## Above the range c is held at its value at the top.
  capped <- day(daily_integration = "one-point", one_point_c_assimilation = 0.7,
                one_point_c_assimilation_slope = 0.15, one_point_c_light_max = 300)
  expect_equal(capped$profit / base$profit, 0.7 + 0.15 * log(300 / 500),
               tolerance = 1e-12)
  ## A zero slope is the constant c, exactly.
  expect_identical(day(daily_integration = "one-point",
                       one_point_c_assimilation = 0.7)$profit,
                   day(daily_integration = "one-point",
                       one_point_c_assimilation = 0.7,
                       one_point_c_light_max = 300)$profit)
  expect_error(day(daily_integration = "one-point",
                   one_point_c_assimilation_slope = 0.1,
                   one_point_c_light_min = 0), "light_min")
  expect_error(day(daily_integration = "one-point",
                   one_point_c_assimilation_slope = 0.1,
                   shading_model = "deep-crown"), "deep-crown")
})

test_that("five-point Lobatto lands on the diurnal integral", {
  ## Three solves; the dark sunrise and sunset nodes are -R_d and no water.
  for (light in c(1, 0.2)) {
    ref <- tf24_day("diurnal", n = 96L, light = light)
    lob <- tf24_day("lobatto", light = light)
    expect_lt(abs(lob$profit / ref$profit - 1), 0.02)
    expect_lt(abs(sum(lob$consumption) / sum(ref$consumption) - 1), 0.02)
  }
})
