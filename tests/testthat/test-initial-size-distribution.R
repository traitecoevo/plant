context("Initial size distribution")

test_that("Multi-species case", {
  p0 <- scm_base_parameters("K93")
  p0$max_patch_lifetime <- 35.10667
  
  # Three species from paper
  sp <- trait_matrix(c(0.042, 0.063, 0.052,
                       8.5e-3, 0.014, 0.015,
                       2.2e-4, 4.6e-4, 3e-4,
                       0.008, 0.008, 0.008,
                       1.8e-4, 4.4e-4, 5.1e-4,
                       1.4e-4, 2.5e-3, 8.8e-3,
                       0.044, 0.044, 0.044),
                     c("b_0", "b_1", "b_2",
                       "c_0", "c_1", "d_0", "d_1"))
  
  p2 <- expand_parameters(sp, p0, mutant = FALSE)
  p2$birth_rate <- c(20, 20, 20)
  
  env <- make_environment("K93")
  ctrl <- scm_base_control()
  
  # manually create an SCM and extract the schedule
  types <- extract_RcppR6_template_types(p2, "Parameters")
  scm <- do.call('SCM', types)(p2, env, ctrl)
  times <- scm$cohort_schedule$all_times
  
  # Create some arbitrary state - while concise, using mapply 
  # fails if not introducing a cohort for one species
  vars <- sapply(scm$state$species, rownames, simplify = F)
  
  state <- mapply(function(names, values, n) 
    matrix(c(values, rep(0, length(names) - 1)),
           nrow=length(names), ncol=n, dimnames=list(names)),
    vars, c(2, 3, 4), c(2, 1, 1), SIMPLIFY=FALSE)

  expect_equal(sapply(state, dim)[2, ], c(2, 1, 1))
  
  # Create introduction schedule for new cohorts
  cohorts <- sapply(state, ncol)
  
  start_time <- sapply(times, function(t) min(t[-1]))
  new_times <- mapply(function(i, t) c(rep(0, i), t[-1]), 
                      cohorts, times, SIMPLIFY = F)
  
  # Run in cohort schedule, then initialse state of initial cohorts
  # this introduces one more cohort than necessary, but if we
  # overwrite the oldest cohorts first then we can just start at t1
  scm$set_cohort_schedule_times(new_times)
  
  # run_next returns species indexes - one cohort for 3 & 2, two cohorts for 1
  expect_equal(scm$run_next(), c(3, 2, 1, 1))
  
  scm$set_state(min(start_time), unlist(z$species), n = cohorts)
  
  # state above + one default cohort (min 2cm height)
  expect_equal(sapply(scm$state$species, rowSums)["height", ], c(6, 5, 6))
  
  # these steps are bundled into make_scm
  scm2 <- make_scm(p2, env, ctrl, state)
  expect_equal(sapply(scm2$state$species, rowSums)["height", ], c(6, 5, 6))
  expect_equal(scm2$cohort_schedule, scm$cohort_schedule)
  
  # verify that this works through the run_scm api
  x <- run_scm_collect(p2, env, ctrl, state)
  
  # check fitness
  expect_equal(x$net_reproduction_ratios, 
               c(3.981e-06, 4.585e-04, 3.514e-04), 
               tolerance = 0.0001)
  
  # check # states, time steps, cohorts
  expect_equal(sapply(x$species, dim)[3, ], c(107, 106, 106))
  
  # check warnings
  expect_error(make_scm(p2, env, ctrl, state[[1]]), 
               "State object has more species than strategies defined in Parameters")
  
  
  # Plot the results
  # t2 <- x$time
  # h1 <- x$species[[1]]["height", , ]
  # h2 <- x$species[[2]]["height", , ]
  # h3 <- x$species[[3]]["height", , ]
  # 
  # cols <- c("#e34a33", "#045a8d", "#000000")
  # 
  # # Species 1 - red
  # matplot(t2, h1, lty=1, col=make_transparent(cols[[1]], .25), type="l",
  #         las=1, xlab="Time (years)", ylab="Height (m)", ylim = c(0, 10))
  # 
  # # Species 2 - blue
  # matlines(t2, h2, lty=1, col=make_transparent(cols[[2]], .25))
  # 
  # # Species 3 - black
  # matlines(t2, h3, lty=1, col=make_transparent(cols[[3]], .25))
})

# Not actually recommending this as a test but it's a good demo for a hard
# problem - starting with just a single cohort makes it difficult to resolve
# whether cohorts should be introduced from the initial size distribution
# or from the cohort schedule.

# `run_scm_error` is somewhat specialised for the build_schedule function,
# particularly persisting the schedule with cohorts at t0. 

# Left unattended, cohorts introduced below the minimum viable height
# of the FF16 strategy trip up the adaptive interpolator. Some ugly patches
# pass the minimum back to split_densities. split_densities doesn't handle
# density properly.
test_that("Build schedule works", {
  env <- make_environment("FF16")
  ctrl <- scm_base_control()
  p0 <- scm_base_parameters("FF16")
  p0$patch_type = 'fixed'
  ctrl$schedule_nsteps = 10
  
  p1 <- expand_parameters(trait_matrix(0.0825, "lma"), p0, FF16_hyperpar,FALSE)
  p1$cohort_schedule_times[[1]] <- seq(0, p1$max_patch_lifetime, length = 10)
  p1$birth_rate <- 20
  
  # manual construction (see also: scm_state)
  init <- matrix(c(10, 0, 0, 0, 0, 0, -3,
                   0.5, 0, 0, 0, 0, 0, -4), ncol = 2)
  rownames(init) <- c("height", "mortality", "fecundity", "area_heartwood",
                      "mass_heartwood", "offspring_produced_survival_weighted",
                      "log_density")

  # interpolate initial size distribution
  splines <- init_spline(list(init), size_idx = 1)
    
  res <- build_schedule(p1, env, ctrl, splines, n_init = 20)
  
  
  expect_false(res$complete)
  expect_equal(res$n_steps, 5)
  expect_equal(attr(res$parameters, "net_reproduction_ratios"), 
               4997.393, tolerance = 1e-6)
  
  # Plot for fun
  par(mfrow=c(2, 1))
  x <- run_scm_collect(res$parameters, env, ctrl)

  t2 <- x$time
  h1 <- x$species[[1]]["height", , ]

  matplot(t2, h1, lty=1, type="l", col = "black",
          las=1, xlab="Time (years)", ylab="Height (m)")
  title("Unresolved schedule")

  x <- run_scm_collect(res$parameters, env, ctrl, res$state)

  t2 <- x$time
  h1 <- x$species[[1]]["height", , ]

  # Species 1
  matplot(t2, h1, lty=1, type="l", col = "black",
          las=1, xlab="Time (years)", ylab="Height (m)")
  title("Resolved schedule")
})
