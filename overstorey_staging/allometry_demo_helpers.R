## Shared computation for TF24_flexible_allometry_demo.qmd (#516).
##
## Kept out of the qmd so it can be exercised by a test rather than only by
## rendering: a chunk that runs is not a chunk that is checked, and prose beside
## a number survives the number being falsified.

## A TF24 strategy with the plasticity switch at `a_pl0` (0 = fixed allometry).
tf24_allometry_strategy <- function(a_pl0 = 1.0, ...) {
  s <- TF24_Strategy(...)
  s$pars$a_pl0 <- a_pl0
  s
}

## A light, well-lit environment with a prescribed soil water content, so the
## demo can drive the carbon balance directly rather than through rainfall and
## the soil's own dynamics.
tf24_demo_env <- function(theta_soil, n_depths = 5L, canopy_top = 40) {
  env <- Environment("TF24")
  env$set_soil_number_of_depths(n_depths)
  env$set_soil_water_state(rep(theta_soil, n_depths))
  env$set_fixed_environment(1.0, canopy_top)
  env
}

## The two departure rates over a range of reserve fractions, at a fixed size.
## Storage is set directly, so `r` is swept rather than emerging -- which is
## what isolates the gate from everything else that moves with it.
allometry_gate_map <- function(a_pl0 = 1.0, height = 5, theta_soil = 0.2,
                               n = 61L) {
  s <- tf24_allometry_strategy(a_pl0)
  env <- tf24_demo_env(theta_soil)

  ## Capacity is needed to turn a target r into a storage value; read it off a
  ## plant seeded at a known fraction of it rather than re-deriving it in R.
  ref <- TF24_Individual(s)
  ref$set_state("height", height)
  ref$set_initial_states(env)
  capacity <- ref$state("storage") / s$pars$a_st3

  r_grid <- seq(0, 1, length.out = n)
  out <- lapply(r_grid, function(r) {
    ind <- TF24_Individual(s)
    ind$set_state("height", height)
    ind$set_state("storage", r * capacity)
    ind$compute_rates(env)
    data.frame(r = r,
               dphi = ind$rate("log_area_leaf_departure"),
               dpsi = ind$rate("log_area_sapwood_departure"),
               dheight = ind$rate("height"))
  })
  cbind(do.call(rbind, out), a_pl0 = a_pl0)
}

## Step one plant through a prescribed soil-water history, so a drought and its
## recovery can both be seen. The ODE is integrated within each interval at
## fixed soil water and the state carried across, which is a piecewise-constant
## driver rather than a smooth one -- adequate for a demo, and stated so the
## reader does not mistake it for how the SCM drives soil water.
allometry_trajectory <- function(a_pl0, soil_history, dt = 0.1, height0 = 4) {
  s <- tf24_allometry_strategy(a_pl0, collect_all_auxiliary = TRUE)
  ind <- TF24_Individual(s)
  ind$set_state("height", height0)
  ind$set_initial_states(tf24_demo_env(soil_history$theta[[1]]))

  nms <- ind$ode_names
  rows <- vector("list", nrow(soil_history))
  t_now <- 0

  for (i in seq_len(nrow(soil_history))) {
    env <- tf24_demo_env(soil_history$theta[[i]])
    span <- soil_history$duration[[i]]
    steps <- max(1L, as.integer(round(span / dt)))
    for (k in seq_len(steps)) {
      y <- ind$ode_state
      res <- grow_individual_to_time(ind, dt, env)
      ind$ode_state <- as.numeric(res$state[1, ])
      t_now <- t_now + dt
    }
    ind$compute_rates(env)
    st <- setNames(as.numeric(ind$ode_state), nms)
    rows[[i]] <- data.frame(
      time = t_now,
      theta_soil = soil_history$theta[[i]],
      height = st[["height"]],
      phi = st[["log_area_leaf_departure"]],
      psi = st[["log_area_sapwood_departure"]],
      survivorship = exp(-st[["mortality"]]),
      area_leaf = ind$aux("competition_effect"),
      area_sapwood = ind$aux("area_sapwood"),
      P = ind$aux("net_mass_production_dt"),
      a_pl0 = a_pl0)
  }
  d <- do.call(rbind, rows)
  ## Two derived quantities the narrative turns on: the canopy relative to what
  ## its height prefers, and the Huber value relative to its preferred one.
  d$canopy_fraction <- exp(d$phi)
  d$huber_ratio <- exp(d$psi)
  d$P_per_area <- d$P / d$area_leaf
  d
}

## A drought with a recovery either side of it.
demo_soil_history <- function(wet = 0.30, dry = 0.11,
                              years_wet = 3, years_dry = 4,
                              step = 0.25) {
  seg <- function(theta, years) {
    n <- as.integer(round(years / step))
    data.frame(theta = rep(theta, n), duration = rep(step, n))
  }
  rbind(seg(wet, years_wet), seg(dry, years_dry), seg(wet, years_wet))
}

## The allometric curves themselves, taken from the strategy's own C++ expansion
## rather than re-derived in R, so a reference line in a plot cannot drift from
## the model it is a reference for.
allometry_reference_curve <- function(s, heights) {
  expand <- get("TF24_strategy_expand_allometry", envir = asNamespace("plant"))
  zero <- rep(0, length(heights))
  z <- expand(s, heights, zero, zero)
  data.frame(height = heights,
             area_leaf = z$area_leaf,
             area_sapwood = z$area_sapwood)
}
