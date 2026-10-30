## Shared computation for TF24_flexible_allometry_demo.qmd (#516).
##
## Kept out of the qmd so it can be exercised by a test rather than only by
## rendering: a chunk that runs is not a chunk that is checked, and prose beside
## a number survives the number being falsified.

## Seat a TF24 individual at a given size and allometric departure, with reserves
## at a_st3 of ITS capacity and the tracked marginal leaf return at its own
## current value -- what a newborn gets at birth size, moved to this size.
##
## set_initial_states() alone cannot do this: it seats a NEWBORN, writing the
## seed's height over whatever height was set before it, so a plant set to 5 m
## and then seeded comes back at birth size.
##
## A copy of tests/testthat/helper-tf24-seat.R, since this file renders outside
## the test suite.
tf24_seat <- function(ind, env, height, phi = 0, psi = 0) {
  ind$set_initial_states(env)
  ind$set_state("height", height)
  ind$set_state("log_area_leaf_departure", phi)
  ind$set_state("log_area_sapwood_departure", psi)
  s <- ind$strategy
  ## mass_sapwood is linear in sapwood area, A*(h) exp(phi) theta exp(psi).
  mass_sapwood <- plant:::TF24_strategy_expand_allometry(s, height, 0, 0)$mass_sapwood *
    exp(phi) * exp(psi)
  ind$set_state("storage", s$pars$a_st3 * s$pars$a_st1 * mass_sapwood)
  ind$compute_rates(env)
  ind$set_state("leaf_marginal_return_tracked", ind$aux("leaf_marginal_return"))
  invisible(ind)
}

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

## The gate's response, swept over SOIL WATER rather than over the reserve pool.
## The gate reads the marginal leaf's carbon balance now (#516), so the pool is
## no longer the axis it responds to -- and holding reserves full while sweeping
## soil water is exactly the demonstration that it is not.
allometry_gate_map <- function(a_pl0 = 1.0, height = 5,
                               theta = seq(0.08, 0.35, length.out = 40L)) {
  s <- tf24_allometry_strategy(a_pl0, collect_all_auxiliary = TRUE)
  out <- lapply(theta, function(th) {
    env <- tf24_demo_env(th)
    ind <- TF24_Individual(s)
    tf24_seat(ind, env, height)      # reserves at a_st3 of capacity throughout
    ind$compute_rates(env)
    data.frame(theta_soil = th,
               margin = ind$aux("leaf_marginal_return"),
               dphi = ind$rate("log_area_leaf_departure"),
               dpsi = ind$rate("log_area_sapwood_departure"),
               dheight = ind$rate("height"))
  })
  cbind(do.call(rbind, out), a_pl0 = a_pl0)
}

## A wet stand left to self-thin, so suppressed plants can be compared with
## dominants. Shade is as much the point as drought: a plant losing the light
## race has the same problem as one losing water, and the same lever.
##
## `amplitude` makes the rain seasonal, rain * (1 + amplitude * sin(2 pi t)),
## floored at zero; `a_sw` switches sapwood acclimation on as well.
allometry_self_thinning <- function(a_pl0 = 1.0, mpl = 25, lma = 0.0825,
                                    birth_rate = 20, rain = 1.5,
                                    amplitude = 0, a_sw = 0) {
  p <- scm_base_parameters("TF24")
  p$max_patch_lifetime <- mpl
  p <- add_strategies(p, trait_matrix(lma, "lma"), hyperpar = TF24_hyperpar,
                      birth_rate = birth_rate)
  p$strategies[[1]]$pars$a_pl0 <- a_pl0
  p$strategies[[1]]$pars$a_sw <- a_sw
  env <- Environment("TF24")
  env$set_soil_number_of_depths(5)
  env$set_soil_water_state(rep(0.428 * 0.5, 5))
  if (amplitude == 0) {
    env$extrinsic_drivers_set_constant("rainfall", rain)
  } else {
    t <- seq(0, mpl, length.out = as.integer(ceiling(mpl * 48)) + 1L)
    env$extrinsic_drivers_set_variable(
      "rainfall", t, pmax(rain * (1 + amplitude * sin(2 * pi * t)), 0))
  }
  out <- run_scm(p, env, Control(), collect = TRUE)
  d <- out$species
  d <- d[is.finite(d$height) & d$density > 0, ]
  d$canopy_fraction <- exp(d$log_area_leaf_departure)
  d$huber_ratio <- exp(d$log_area_sapwood_departure)
  d$survivorship <- exp(-d$mortality)
  d$a_pl0 <- a_pl0
  list(cohorts = d, R0 = out$offspring_production)
}

## Step one plant through a prescribed soil-water history, so a drought and its
## recovery can both be seen. The ODE is integrated within each interval at
## fixed soil water and the state carried across, which is a piecewise-constant
## driver rather than a smooth one -- adequate for a demo, and stated so the
## reader does not mistake it for how the SCM drives soil water.
allometry_trajectory <- function(a_pl0, soil_history, dt = 0.1, height0 = 4) {
  s <- tf24_allometry_strategy(a_pl0, collect_all_auxiliary = TRUE)
  ind <- TF24_Individual(s)
  tf24_seat(ind, tf24_demo_env(soil_history$theta[[1]]), height0)

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

## The carbon budget split by whether a cost scales with leaf area, which is what
## the shedding criterion rests on. Returns the three terms plus the elasticity
## of per-leaf assimilation to hydraulic supply.
##
## `eta` is estimated by differencing at FIXED sapwood area -- lower the leaf
## departure by d and raise the sapwood departure by d, so A_s is unchanged. That
## is what real thinning does; letting A_s follow A down the pipe model measures
## something else entirely.
allometry_carbon_terms <- function(height, theta_soil, d = 0.02, phi = 0,
                                   psi = 0) {
  p <- TF24_Strategy()$pars
  eta_c <- 1 - 2 / (1 + p$eta) + 1 / (1 + 2 * p$eta)
  cc <- p$a_bio * p$a_y
  CONV <- 60 * 60 * 12 * 365 / 1e6          # per-second to annual, mol
  s <- TF24_Strategy(collect_all_auxiliary = TRUE)

  at <- function(ph, ps) {
    e <- tf24_demo_env(theta_soil)
    i <- TF24_Individual(s)
    tf24_seat(i, e, height, phi = ph, psi = ps)
    i$compute_rates(e)
    A <- i$aux("competition_effect")
    list(A = A, A_s = i$aux("area_sapwood"),
         abar = (i$aux("profit") + i$aux("shadow_cost")) * CONV,
         P = i$aux("net_mass_production_dt"))
  }

  a <- at(phi, psi)
  b <- at(phi - d, psi + d)                 # same A_s, less A
  stopifnot(abs(b$A_s / a$A_s - 1) < 1e-9)

  ## Costs that scale with leaf area: leaf, fine root, and bark (pinned to leaf
  ## area). Sapwood is the one that does not.
  per_area <- cc * (p$r_l * p$lma + p$r_r * p$a_r1 +
                      p$r_b * p$a_b1 * p$theta * height * eta_c * p$rho) +
    (p$k_l * p$lma + p$k_r * p$a_r1 +
       p$k_b * p$a_b1 * p$theta * height * eta_c * p$rho)
  m_s <- a$A_s * height * eta_c * p$rho
  dPdA <- (b$P - a$P) / (b$A - a$A)

  ## eta is recovered from the measured dP/dA and the (exact) decomposition,
  ## NOT by differencing abar. The profit auxes misreport assimilation off the
  ## trajectory by ~2 per cent, which is small against the level but comparable
  ## to the difference being taken, so a direct estimate of eta is unusable
  ## while dP/dA -- a difference of P itself -- is not.
  list(height = height, theta_soil = theta_soil, psi = psi,
       abar = a$abar,
       gain = cc * a$abar * a$A,
       g = per_area * a$A,
       s = cc * p$r_s * m_s + p$k_s * m_s,
       P = a$P,
       kappa = per_area / (cc * a$abar),
       eta = 1 - (dPdA + per_area) / (cc * a$abar),
       dPdA = dPdA)
}

## The sapwood controller's signal, swept over the Huber value at fixed height
## and leaf area (#516). Returns the marginal growth return R_s beside the growth
## rate it is the derivative of, so the two can be plotted on the same axis and
## the zero-crossing checked against the peak rather than asserted.
allometry_sapwood_controller <- function(theta_soil, height = 10, a_sw = 1.0,
                                         psi = seq(-0.3, 0.7, length.out = 41L)) {
  s <- TF24_Strategy(collect_all_auxiliary = TRUE)
  s$pars$a_sw <- a_sw
  out <- lapply(psi, function(p) {
    env <- tf24_demo_env(theta_soil)
    ind <- TF24_Individual(s)
    tf24_seat(ind, env, height, psi = p)
    ind$compute_rates(env)
    data.frame(psi = p, huber_ratio = exp(p),
               R_s = ind$aux("sapwood_marginal_return"),
               dheight = ind$rate("height"),
               P = ind$aux("net_mass_production_dt"))
  })
  d <- do.call(rbind, out)
  d$theta_soil <- theta_soil
  d$height <- height
  d
}

## Step one plant: acclimate wet, zero the mortality integral, drought, recover.
## Zeroing mortality is what makes sizes comparable -- otherwise a tall plant's
## answer is dominated by the lifetime of hazard it took to get tall, and the
## drought's own contribution is unreadable.
allometry_drought <- function(a_pl0, height, theta_dry, years_dry = 4,
                              wet = 0.30, acclimate = 2, recover = 4,
                              a_sw = 0, dt = 0.1) {
  s <- TF24_Strategy(collect_all_auxiliary = TRUE)
  s$pars$a_pl0 <- a_pl0
  s$pars$a_sw  <- a_sw
  ind <- TF24_Individual(s)
  tf24_seat(ind, tf24_demo_env(wet), height)
  ## The canopy's lowest point, read at every step: by the end of the recovery
  ## a surviving plant has rebuilt, so the final state would hide the thinning.
  min_canopy <- 1
  run <- function(theta, years) {
    env <- tf24_demo_env(theta)
    for (k in seq_len(as.integer(round(years / dt)))) {
      res <- try(grow_individual_to_time(ind, dt, env), silent = TRUE)
      if (inherits(res, "try-error")) return(FALSE)
      ind$ode_state <- as.numeric(res$state[1, ])
      min_canopy <<- min(min_canopy,
                         exp(ind$state("log_area_leaf_departure")))
    }
    TRUE
  }
  if (!run(wet, acclimate)) return(NULL)
  ind$set_state("mortality", 0)
  if (!run(theta_dry, years_dry)) return(NULL)
  if (!run(wet, recover)) return(NULL)
  st <- setNames(as.numeric(ind$ode_state), ind$ode_names)
  data.frame(a_pl0 = a_pl0, height = height, theta_dry = theta_dry,
             years_dry = years_dry,
             survivorship = exp(-st[["mortality"]]),
             canopy = exp(st[["log_area_leaf_departure"]]),
             min_canopy = min_canopy,
             huber = exp(st[["log_area_sapwood_departure"]]))
}

## The intensity-by-size grid the "does not buy survival" section rests on.
allometry_drought_grid <- function(heights = c(10, 20, 24),
                                   thetas = c(0.145, 0.140, 0.135),
                                   years_dry = 4) {
  g <- expand.grid(height = heights, theta_dry = thetas)
  out <- lapply(seq_len(nrow(g)), function(i) {
    a <- allometry_drought(0.0, g$height[i], g$theta_dry[i], years_dry)
    b <- allometry_drought(1.0, g$height[i], g$theta_dry[i], years_dry)
    if (is.null(a) || is.null(b)) return(NULL)
    data.frame(height = g$height[i], theta_dry = g$theta_dry[i],
               S_fixed = a$survivorship, S_shed = b$survivorship,
               min_canopy = b$min_canopy,
               gain = b$survivorship / a$survivorship)
  })
  do.call(rbind, out)
}

## Maintenance -- respiration in carbon units plus turnover, fully charged -- by
## tissue, as shares of the whole bill, for a plant on its allometry. Which
## shares follow the canopy (leaf, fine root, bark) and which do not (sapwood) is
## what bounds what thinning can save.
allometry_cost_shares <- function(height) {
  s <- TF24_Strategy()
  p <- s$pars
  z <- plant:::TF24_strategy_expand_allometry(s, height, 0, 0)
  cc <- p$a_bio * p$a_y
  bill <- c(leaf = cc * p$r_l * z$mass_leaf + p$k_l * z$mass_leaf,
            bark = cc * p$r_b * z$mass_bark + p$k_b * z$mass_bark,
            sapwood = cc * p$r_s * z$mass_sapwood + p$k_s * z$mass_sapwood,
            root = cc * p$r_r * z$mass_root + p$k_r * z$mass_root)
  bill / sum(bill)
}

## Net production per unit leaf area as a plant thins at fixed height, with the
## stem kept at 0.35 of the leaf area it loses -- sapwood lost to heartwood
## while the canopy goes faster -- so stem per leaf rises as it would.
allometry_thinning_per_area <- function(height = 10, theta_soil = 0.13,
                                        phi = c(0, -0.7, -1.6, -2.3, -3.0,
                                                -4.0, -5.0)) {
  out <- lapply(phi, function(k) {
    ind <- TF24_Individual(TF24_Strategy(collect_all_auxiliary = TRUE))
    env <- tf24_demo_env(theta_soil)
    tf24_seat(ind, env, height, phi = k, psi = -k * 0.35)
    ind$compute_rates(env)
    data.frame(phi = k, canopy = exp(k), huber_ratio = exp(-k * 0.35),
               P_per_area = ind$aux("net_mass_production_dt") /
                 ind$aux("competition_effect"))
  })
  do.call(rbind, out)
}

## One shedding plant through a drought, step by step: acclimate wet, zero the
## mortality integral, then record the canopy and survivorship every `dt` of the
## drought. What it shows is the ORDER of the two -- whether a plant thins before
## or after its survivorship collapses -- which an end-of-run summary cannot.
allometry_drought_course <- function(height = 10, theta_dry = 0.140, years = 2,
                                     wet = 0.30, acclimate = 2, dt = 0.1) {
  ind <- TF24_Individual(tf24_allometry_strategy(1.0, collect_all_auxiliary = TRUE))
  tf24_seat(ind, tf24_demo_env(wet), height)
  step <- function(env) {
    r <- grow_individual_to_time(ind, dt, env)
    ind$ode_state <- as.numeric(r$state[1, ])
  }
  for (k in seq_len(round(acclimate / dt))) step(tf24_demo_env(wet))
  ind$set_state("mortality", 0)
  dry <- tf24_demo_env(theta_dry)
  out <- lapply(seq_len(round(years / dt)), function(k) {
    step(dry)
    data.frame(year = k * dt,
               canopy = exp(ind$state("log_area_leaf_departure")),
               survivorship = exp(-ind$state("mortality")))
  })
  do.call(rbind, out)
}
