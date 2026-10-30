  ## Three years of seasonal rain, on most of whose steps the young stand takes
  ## up a small fraction of the water flowing through the soil. Each step keeps
  ## its record, so a pinned run, the sweep and an invader take it as the run did.
  p0 <- scm_base_parameters("TF24", "TF24_Env")
  p0$max_patch_lifetime <- 3
  p0$patch_type <- "fixed"
  lma <- 0.1978791
  p <- add_strategies(p0, trait_matrix(lma, "lma"))
  env <- function() {
    e <- Environment("TF24")
    t <- seq(0, 3, length.out = 601)
    e$extrinsic_drivers_set_variable("rainfall", t,
                                     0.25 * (1 + 0.95 * sin(2 * pi * t)))
    e
  }
  ctrl <- Control(node_density_in_birth_date = TRUE, ode_tol_abs = 1e-8,
                  ode_soil_substep = TRUE)
  substepped <- run_scm(p, env(), ctrl, keep_states = TRUE)
  expect_gt(mean(lengths(substepped$ode_soil_substep_ends) > 0), 0.5)

  ## A pinned run of the schedule, records and all, repeats the run bit for bit.
  rerun <- function(q) {
    scm <- SCM("TF24", "TF24_Env")(q, env(), empty_events(), ctrl)
    sched <- scm$node_schedule
    sched$all_times <- substepped$node_schedule$all_times
    sched$set_ode_steps(substepped$ode_times, substepped$ode_step_sizes,
                        substepped$ode_soil_input_slopes,
                        substepped$ode_soil_substep_ends)
    scm$node_schedule <- sched
    scm$run()
    scm
  }
  again <- rerun(p)
  expect_identical(again$patch$ode_state, substepped$patch$ode_state)
  expect_identical(again$ode_soil_substep_ends, substepped$ode_soil_substep_ends)
  ## So does the schedule carried in the parameters. Its times alone, a grid,
  ## are refused where a step would substep the soil, and its sizes without
  ## their records where a step spans a knot of the rain.
  q <- p
  q$ode_times <- substepped$ode_times
  expect_error(run_scm(q, env(), ctrl), "step the run's schedule")
  q$ode_step_sizes <- substepped$ode_step_sizes
  expect_error(run_scm(q, env(), ctrl), "carries no substeps")
  q$ode_soil_input_slopes <- substepped$ode_soil_input_slopes
  q$ode_soil_substep_ends <- substepped$ode_soil_substep_ends
  carried <- run_scm(q, env(), ctrl)
  expect_identical(carried$patch$ode_state, substepped$patch$ode_state)

  ## The sweep holds each step's substeps and input slopes, as a tangent of the
  ## pinned schedule does, and both agree with its central differences.
  g <- stand_gradient(substepped, traits = "1.lma", hyperpar = NULL)
  traits <- census_trait_names_tf24(substepped)
  tangent <- ladder_trajectory_tangent_tf24(substepped,
                                            as.numeric(traits == "1.lma"))
  ## The sweep's rows: the census, then each species' fitness.
  metrics <- function(scm) c(stand_census(scm), stand_fitness(scm))
  expect_identical(names(metrics(substepped)), rownames(g$gradient))
  ## The tangent is the run at another scalar, not the run, so its value agrees
  ## to the budget test-gradient-ladder-introductions.R explains, not bit for bit
  ## (3.6e-13 on macOS, 1.25e-12 on Windows).
  expect_equal(tangent$value, unname(metrics(substepped)), tolerance = 1e-9)
  expect_equal(unname(g$gradient[, "1.lma"]), tangent$tangent, tolerance = 1e-6)
  at_lma <- function(x) {
    q <- p
    s <- q$strategies[[1]]
    pars <- s$pars
    pars$lma <- x
    s$pars <- pars
    q$strategies[[1]] <- s
    metrics(rerun(q))
  }
  ## The difference's truncation error at a relative step of 1e-4 is 7% on the
  ## first row; at 1e-5 it is 0.4%, and at 1e-6 rounding takes over.
  d <- 1e-5 * lma
  differenced <- (at_lma(lma + d) - at_lma(lma - d)) / (2 * d)
  expect_equal(unname(g$gradient[, "1.lma"]), unname(differenced),
               tolerance = 1e-2)

  ## An invader with the stand's traits, integrated on the run's steps with the
  ## soil read from the run, has the stand's fitness. Not to the bit: on this
  ## stand it misses by 4e-14 with the soil coupled too.
  ratios <- substepped$net_reproduction_ratios
  substepped$run_mutant(p)
  expect_equal(substepped$net_reproduction_ratios, ratios, tolerance = 1e-12)
