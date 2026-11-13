## Runs that use none of the opt-in features, and what each reports. Shared by
## test-default-baseline.R and scripts/capture-default-baseline.R. The runs go
## through run() (#428), so the capture needs a build that has it; the recorded
## baseline came from a build before #661, through the SCM entry point run()
## replaced, which drove the solver the same way.
default_baseline_cases <- function() {
  outputs <- function(scm) {
    list(ode_state = scm$patch$ode_state,
         ode_times = scm$ode_times,
         net_reproduction_ratios = scm$net_reproduction_ratios,
         offspring_production = scm$offspring_production)
  }
  ff16 <- function(lifetime) {
    p <- testing_parameters("FF16")
    p$max_patch_lifetime <- lifetime
    p$node_schedule_times <- list()
    add_strategies(p, trait_matrix(0.0825, "lma"))
  }
  out <- list()

  ## FF16 and K93 on the defaults: Weibull patch ages, the height coordinate.
  out$ff16 <- outputs(run(ff16(105.32))$solver)
  p <- testing_parameters("K93")
  p <- add_strategies(p, trait_matrix(0.059, "b_0"), birth_rate = 20)
  out$k93 <- outputs(run(p)$solver)

  ## FF16 on the birth-date coordinate with a harvest, no disturbance regime.
  p <- ff16(40)
  p$patch_type <- "fixed"
  ev <- events(node_introductions(p), harvest(time = 20, fraction = 0.3))
  out$ff16_birth_date_harvest <-
    outputs(run(p, events = ev,
                control = Control(node_density_coordinate = "birth_date"))$solver)

  ## An invasion of the default FF16 stand.
  scm <- run(ff16(105.32))$solver
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
  first <- run(fixed(10))$solver
  state <- export_patch_state(first)
  p2 <- set_initial_state(fixed(20), state)
  tt <- fixed(20)$node_schedule_times[[1]]
  p2$node_schedule_times <- list(tt[tt > state$time])
  out$ff16_resumed <- outputs(run(p2)$solver)

  ## A short TF24 stand under constant rain, on the birth-date coordinate.
  p <- testing_parameters("TF24")
  p$max_patch_lifetime <- 5
  p <- add_strategies(p, trait_matrix(0.0825, "lma"), hyperpar = default_hyperpar("TF24"),
                      birth_rate = 1)
  env <- Environment("TF24")
  env$extrinsic_drivers_set_constant("rainfall", 1)
  out$tf24 <- outputs(run(p, env, Control(node_density_coordinate = "birth_date",
                                          ode_step_size_max = 0.05))$solver)
  out
}
