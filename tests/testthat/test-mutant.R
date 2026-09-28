## What an invasion run has to reproduce. Every number here was measured on
## develop, against a recorder reached through three hooks the ODE solver called
## into the patch. odelia's rewrite deleted those hooks; what replaces them is
## odelia's own store/load channel -- the run keeps the field in the same
## per-(step, stage) row it already keeps what a rate evaluation solved for, and
## the invasion pass loads it. The numbers did not move, which is the point of
## keeping them: they were the specification the replacement was written to, not
## a re-pin taken from it.

test_that("mutant method works", {
  # basic setup
  p0 <- scm_base_parameters("FF16")
  p0$max_patch_lifetime <- 50

  e <- Environment("FF16")
  ctrl <- Control()

  tol <- 1e-4
  
  # We'll run tests with 1 and 3 residents, each with different numbers of mutants
  
  lma <- c(0.05, 0.1, 0.2)
  birth_rate <- 1

  # 1 resident strategies
  pr1 <- add_strategies(p0, trait_matrix(lma[2], "lma"), birth_rate = rep(birth_rate, 1))

  pr1m1 <- add_strategies(pr1, trait_matrix(lma[3], "lma"), birth_rate = rep(birth_rate, 1))

  pr1m3 <- add_strategies(pr1, trait_matrix(lma, "lma"), birth_rate = rep(birth_rate, 3))

  pr1m10 <- add_strategies(pr1, trait_matrix(seq(lma[1], lma[3], length.out=10), "lma"), birth_rate = rep(birth_rate, 10))

  # test error handling
  # scm object but not yet run
  types <- extract_RcppR6_template_types(pr1, "Parameters")
  scm <- do.call("SCM", types)(pr1, e, empty_events(), ctrl)

  expect_error(scm$run_mutant(p0), "Run a resident first")

  # check mutant fitness against resindet and expected values
  scm <- run_scm(pr1, e, ctrl)
  pr1_rr <- scm$net_reproduction_ratios
  expected <- 2.77322
  expect_equal(pr1_rr, expected, tolerance = tol)

  scm$run_mutant(pr1m1)
  pr1m1_rr <- scm$net_reproduction_ratios
  expected <- c(2.77322, 3.707605)
  expect_equal(pr1m1_rr, expected, tolerance = tol)
  expect_equal(pr1m1_rr[1], pr1_rr, tolerance = tol)

  scm$run_mutant(pr1m3)
  pr1m3_rr <- scm$net_reproduction_ratios
  expected <- c(2.77322, 3.7429e-10, 2.77322, 3.70753)
  expect_equal(pr1m3_rr, expected, tolerance = tol)
  expect_equal(pr1m3_rr[1], pr1_rr, tolerance = tol)

  scm$run_mutant(pr1m10)
  pr1m10_rr <- scm$net_reproduction_ratios
  expected <- c(2.773222, 3.742935e-10, 9.308944e-07, 0.1363641, 2.773222, 3.890554, 1.524582, 1.160212, 1.871261, 2.765328, 3.707372)
  expect_equal(pr1m10_rr, expected, tolerance = tol)
  expect_equal(pr1m10_rr[1], pr1_rr, tolerance = tol)

  # 3 resident strategies
  pr3 <- add_strategies(p0, trait_matrix(lma, "lma"), birth_rate = rep(birth_rate, 3))
  
  pr3m1 <- add_strategies(pr3, trait_matrix(lma[3], "lma"), birth_rate = rep(birth_rate, 1))

  pr3m3 <- add_strategies(pr3, trait_matrix(lma, "lma"), birth_rate = rep(birth_rate, 3))

  pr3m10 <- add_strategies(pr3, trait_matrix(seq(lma[1], lma[3], length.out = 10), "lma"), birth_rate = rep(birth_rate, 10))

  scm <- run_scm(pr3, e, ctrl)
  pr3_rr <- scm$net_reproduction_ratios
  expected <- c(4.265e-10, 2.831741, 0.09125339)
  expect_equal(pr3_rr, expected, tolerance = tol)


  scm$run_mutant(pr3m1)
  pr3m1_rr <- scm$net_reproduction_ratios
  expected <- c(4.265e-10, 2.831741, 0.09125339, 0.09125339)
  expect_equal(pr3m1_rr, expected, tolerance = tol)
  expect_equal(pr3m1_rr[1:3], pr3_rr, tolerance = tol)

  scm$run_mutant(pr3m3)
  pr3m3_rr <- scm$net_reproduction_ratios
  expected <- c(4.265e-10, 2.831741, 0.09125339, 4.265e-10, 2.831741, 0.09125339)
  expect_equal(pr3m3_rr, expected, tolerance = tol)
  expect_equal(pr3m3_rr[1:3], pr3_rr, tolerance = tol)

  scm$run_mutant(pr3m10)
  pr3m10_rr <- scm$net_reproduction_ratios
  expected <- c(4.265011e-10, 2.831741, 0.09125377, 4.265011e-10, 5.587752e-06, 0.266188, 2.831741, 2.690585, 0.3796333, 0.07098642, 0.07226859, 0.08342181, 0.09125377)
  expect_equal(pr3m10_rr, expected, tolerance = tol)
  expect_equal(pr3m3_rr[1:3], pr3_rr, tolerance = tol)
})

test_that("mutant method densities", {
  # For a mutant strategy identical to the resident, the mutant method must
  # reproduce exactly the fitness that strategy attains when run as a resident.
  # This is an identity of the machinery rather than a near-equilibrium
  # approximation -- it holds at any birth rate and any patch lifetime (the two
  # sides agree to ~1e-13 below, far inside the 1e-3 tolerance). We therefore
  # check the invariant across a spread of birth rates (including the degenerate
  # zero-birth case) at two patch lifetimes.
  #
  # Short patch lifetimes are used on purpose: because the agreement is
  # lifetime-independent, a shorter patch retains the full strength of the check
  # while running several times faster than the model's default lifetime. (The
  # earlier versions sampled birth rates around a hard-coded equilibrium at the
  # default lifetime, but the test never asserted anything *about* that
  # equilibrium -- only the resident-vs-mutant identity -- so the long, costly
  # patch bought no extra coverage.)
  ctrl <- Control()

  traits <- trait_matrix(0.0825, c("lma"))
  tol <- 1e-3

  # fitness at birth rate x computed two ways: as a resident, and as a mutant
  # of the resident -- which must agree. The identity holds against whatever
  # competitive landscape the resident run produces, so we do NOT refine the
  # cohort schedule first: refinement is irrelevant to the invariant but was the
  # dominant cost (it roughly tripled this block's run time).
  f_test <- function(p, x) {
    p1 <- p
    p1$strategies[[1]]$birth_rate_y <- x

    scm <- run_scm(p1, ctrl = ctrl)
    r_rr <- scm$net_reproduction_ratios

    scm$run_mutant(p1)
    m_rr <- scm$net_reproduction_ratios

    dplyr::tibble(birth_rate = x, resident_f = log(r_rr), mutant_f = log(m_rr))
  }

  run_case <- function(life, birth_rates) {
    p0 <- scm_base_parameters("FF16")
    p0$max_patch_lifetime <- life
    pr1 <- add_strategies(p0, traits, birth_rate = 1)

    outputs <- purrr::map_df(birth_rates, ~ f_test(pr1, .x))

    expect_equal(birth_rates, outputs$birth_rate, tolerance = tol)
    expect_equal(outputs$resident_f, outputs$mutant_f, tolerance = tol)
  }

  run_case(30, c(0, 5, 10, 20))
  run_case(20, c(0, 5, 10, 20))
})

test_that("mutant method densities, TF24", {
  # The same resident-vs-mutant identity as the block above, for a model whose
  # rates refuse a state. TF24's storage pool reports an overshoot below empty by
  # throwing a domain error, which the adaptive stepper answers by shrinking and
  # retrying -- a routine event, some hundreds of times in a resident run that
  # goes on to complete normally.
  #
  # run_mutant() pins the stepper to the resident's recorded step times, and that
  # path used to call the stepper with no domain handling at all, so the first of
  # those refusals killed the replay (#642). It failed for every TF24 strategy
  # tried, this identity case included.
  #
  # ⚠️ THE FIX IS NOT THE ONE THIS COMMENT USED TO NAME. odelia 0.4.0 subdivided a
  # refused pinned step; this replay does not subdivide at all, deliberately --
  # see NEWS.md under Known issues, where the reason is that subdividing is what
  # made the answer depend on how many invaders shared the call. A replay that met
  # a refusal would now FAIL rather than shrink, and `expect_no_error` below
  # passing is the statement that it meets none.
  #
  # ⚠️ THE LIFETIME BUYS THE REGIME AND THE COHORT COUNT IS WHAT COSTS, and the
  # default schedule confounds them by deriving its introduction count from the
  # lifetime. Held at 6 and varying only how many of that schedule's 89
  # introductions are kept, against the 1e-3 this compares at:
  #
  #   introductions     20     40     60     89
  #   steps            311    363    452    497
  #   log gap        4e-15  4e-13  1e-14  4e-15
  #   seconds          2.7    6.0   11.2     18
  #
  # The gap does not fall with either count, because this is an identity rather
  # than an approximation. So a longer recording buys only more of the regime,
  # and twenty introductions are enough to hold the identity at 4e-15.
  #
  # ⚠️ THOSE STEP COUNTS ARE A TWENTIETH OF WHAT THEY WERE, and the pool is why.
  # This table read 6071 to 12714 steps and 71 to 623 seconds while compute_rates
  # ran the clamped pre-v9 pool the templating commit transcribed; restoring the
  # charge and drain form took the same four fixtures to 311 to 497. The stiffness
  # the schedule was thinned to avoid was mostly the pool integrating past its own
  # ceiling. Both counts are asserted below rather than left to the constants.
  ctrl <- Control()
  tol <- 1e-3

  p0 <- scm_base_parameters("TF24")
  p0$max_patch_lifetime <- 6
  p1 <- add_strategies(p0, trait_matrix(0, "TF24_floor_lambda_o"),
                       hyperpar = TF24_hyperpar, birth_rate = 1)
  full <- p1$node_schedule_times[[1]]
  p1$node_schedule_times <-
    list(full[round(seq(1, length(full), length.out = 20))])

  env <- Environment("TF24")
  env$set_soil_water_state(rep(0.428 * 0.5, env$get_soil_number_of_depths()))
  env$extrinsic_drivers_set_constant("rainfall", 1)

  scm <- run_scm(p1, env = env, ctrl = ctrl)
  resident_rr <- scm$net_reproduction_ratios

  # Not assertions about the model, just guards that the replay below has
  # something to replay: a resident that died out would make the identity
  # trivial, and a short recording would make it cheap in the wrong way.
  expect_true(all(is.finite(resident_rr)) && all(resident_rr > 0))
  expect_gt(length(scm$ode_times), 250)
  expect_equal(scm$patch$species[[1]]$size, 20L)

  # Identical mutant, replaying the resident's own recorded environment, must
  # recover the resident's own fitness.
  expect_no_error(scm$run_mutant(p1))
  mutant_rr <- scm$net_reproduction_ratios

  expect_equal(log(mutant_rr), log(resident_rr), tolerance = tol)
})
