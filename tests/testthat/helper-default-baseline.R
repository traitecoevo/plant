## Runs that use none of the opt-in features, and what each reports. Shared by
## test-default-baseline.R and scripts/capture-default-baseline.R, which
## records them from a build that predates the features. Uses only API that
## predates #661, so the capture can run against such a build.
default_baseline_cases <- function() {
  outputs <- function(scm) {
    list(ode_state = scm$patch$ode_state,
         ode_times = scm$ode_times,
         net_reproduction_ratios = scm$net_reproduction_ratios,
         offspring_production = scm$offspring_production)
  }
  ff16 <- function(lifetime) {
    p <- scm_base_parameters("FF16")
    p$max_patch_lifetime <- lifetime
    p$node_schedule_times <- list()
    add_strategies(p, trait_matrix(0.0825, "lma"))
  }
  out <- list()

  ## FF16 and K93 on the defaults: Weibull patch ages, the height coordinate.
  out$ff16 <- outputs(run_scm(ff16(105.32)))
  p <- scm_base_parameters("K93")
  p <- add_strategies(p, trait_matrix(0.059, "b_0"), birth_rate = 20)
  out$k93 <- outputs(run_scm(p))

  ## FF16 on the birth-date coordinate with a harvest, no disturbance regime.
  p <- ff16(40)
  p$patch_type <- "fixed"
  ev <- events(node_introductions(p), harvest(time = 20, fraction = 0.3))
  out$ff16_birth_date_harvest <-
    outputs(run_scm(p, events = ev, ctrl = Control(node_density_coordinate = "birth_date")))

  ## An invasion of the default FF16 stand.
  scm <- run_scm(ff16(105.32))
  scm$run_mutant(add_mutant(ff16(105.32), trait_matrix(0.1, "lma")))
  out$ff16_invader <- outputs(scm)

  ## A resume of an FF16 stand from its export, under one regime throughout: a
  ## resume re-stamps birth weights from the resuming regime (#661), which under
  ## the Weibull moves with the lifetime.
  fixed <- function(lifetime) {
    p <- ff16(lifetime)
    p$patch_type <- "fixed"
    p
  }
  first <- run_scm(fixed(10))
  state <- export_patch_state(first)
  p2 <- set_initial_state(fixed(20), state)
  tt <- fixed(20)$node_schedule_times[[1]]
  p2$node_schedule_times <- list(tt[tt > state$time])
  out$ff16_resumed <- outputs(run_scm(p2))

  ## A short TF24 stand under constant rain, on the birth-date coordinate.
  p <- scm_base_parameters("TF24")
  p$max_patch_lifetime <- 5
  p <- add_strategies(p, trait_matrix(0.0825, "lma"), hyperpar = TF24_hyperpar,
                      birth_rate = 1)
  env <- Environment("TF24")
  env$extrinsic_drivers_set_constant("rainfall", 1)
  out$tf24 <- outputs(run_scm(p, env, Control(node_density_coordinate = "birth_date",
                                              ode_step_size_max = 0.05)))
  out
}
