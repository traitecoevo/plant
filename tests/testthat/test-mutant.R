context("Mutant-method")

test_that("mutant method works", {
  # basic setup 
  p0 <- scm_base_parameters("FF16")
  p0$max_patch_lifetime <- 50
  
  e <- make_environment("FF16")
  ctrl <- scm_base_control()
  ctrl$save_RK45_cache = TRUE
    
  tol <- 1e-4
  
  # We'll run tests with 1 and 3 residents, each with different numbers of mutants
  
  lma <- c(0.05, 0.1, 0.2)
  birth_rate <- 1

  # 1 resident strategies
  pr1 <- expand_parameters(trait_matrix(lma[2], "lma"), p0,
    birth_rate_list = rep(birth_rate, 1)
  )

  pr1m1 <- expand_parameters(trait_matrix(lma[3], "lma"), pr1,
    birth_rate_list = rep(birth_rate, 1)
  )

  pr1m3 <- expand_parameters(trait_matrix(lma, "lma"), pr1,
    birth_rate_list = rep(birth_rate, 3)
  )

  pr1m10 <- expand_parameters(trait_matrix(seq(lma[1], lma[3], length.out=10), "lma"), pr1,
    birth_rate_list = rep(birth_rate, 10)
  )

  # test error handling
  # scm object but not yet run
  types <- extract_RcppR6_template_types(pr1, "Parameters")
  scm <- do.call("SCM", types)(pr1, e, ctrl)

  expect_error(scm$run_mutant(p0), "Run a resident first to generate a competitve landscape") 

  # check mutant fitness against resindet and expected values
  scm <- run_scm(pr1, e, ctrl)
  pr1_rr <- scm$net_reproduction_ratios
  expected <- 2.77322
  expect_equal(pr1_rr, expected, tol = tol)

  scm$run_mutant(pr1m1)
  pr1m1_rr <- scm$net_reproduction_ratios
  expected <- c(2.77322, 3.707605)
  expect_equal(pr1m1_rr, expected, tol = tol)
  expect_equal(pr1m1_rr[1], pr1_rr, tol = tol)

  scm$run_mutant(pr1m3)
  pr1m3_rr <- scm$net_reproduction_ratios
  expected <- c(2.77322, 3.7429e-10, 2.77322, 3.70753)
  expect_equal(pr1m3_rr, expected, tol = tol)
  expect_equal(pr1m3_rr[1], pr1_rr, tol = tol)

  scm$run_mutant(pr1m10)
  pr1m10_rr <- scm$net_reproduction_ratios
  expected <- c(2.773222, 3.742935e-10, 9.308944e-07, 0.1363641, 2.773222, 3.890554, 1.524582, 1.160212, 1.871261, 2.765328, 3.707372)
  expect_equal(pr1m10_rr, expected, tol = tol)
  expect_equal(pr1m10_rr[1], pr1_rr, tol = tol)

  # 3 resident strategies
  pr3 <- expand_parameters(trait_matrix(lma, "lma"), p0,
    birth_rate_list = rep(birth_rate, 3)
  )
  
  pr3m1 <- expand_parameters(trait_matrix(lma[3], "lma"), pr3,
    birth_rate_list = rep(birth_rate, 1)
  )

  pr3m3 <- expand_parameters(trait_matrix(lma, "lma"), pr3,
    birth_rate_list = rep(birth_rate, 3)
  )

  pr3m10 <- expand_parameters(trait_matrix(seq(lma[1], lma[3], length.out = 10), "lma"), pr3,
    birth_rate_list = rep(birth_rate, 10)
  )

  scm <- run_scm(pr3, e, ctrl)
  pr3_rr <- scm$net_reproduction_ratios
  expected <- c(4.265e-10, 2.831741, 0.09125339)
  expect_equal(pr3_rr, expected, tol = tol)


  scm$run_mutant(pr3m1)
  pr3m1_rr <- scm$net_reproduction_ratios
  expected <- c(4.265e-10, 2.831741, 0.09125339, 0.09125339)
  expect_equal(pr3m1_rr, expected, tol = tol)
  expect_equal(pr3m1_rr[1:3], pr3_rr, tol = tol)

  scm$run_mutant(pr3m3)
  pr3m3_rr <- scm$net_reproduction_ratios
  expected <- c(4.265e-10, 2.831741, 0.09125339, 4.265e-10, 2.831741, 0.09125339)
  expect_equal(pr3m3_rr, expected, tol = tol)
  expect_equal(pr3m3_rr[1:3], pr3_rr, tol = tol)

  scm$run_mutant(pr3m10)
  pr3m10_rr <- scm$net_reproduction_ratios
  expected <- c(4.265011e-10, 2.831741, 0.09125377, 4.265011e-10, 5.587752e-06, 0.266188, 2.831741, 2.690585, 0.3796333, 0.07098642, 0.07226859, 0.08342181, 0.09125377)
  expect_equal(pr3m10_rr, expected, tol = tol)
  expect_equal(pr3m3_rr[1:3], pr3_rr, tol = tol)
})

test_that("mutant method densities", {
  # Now test with different resident densities
    
  # 1 resident strategies
  ctrl <- scm_base_control()
  ctrl$save_RK45_cache = TRUE

  lma_attr <- 0.0825
  p0 <- scm_base_parameters("FF16")
  traits <- trait_matrix(lma_attr, c("lma"))
  tol <- 1e-3

  # Function calculates fitness at different brith_rates, using
  # 2 methods: as resident, as mutant
  # These should all agree
  f_test <- function(p, x, traits) {
    p1 <- p
    p1$strategies[[1]]$birth_rate_y <- x

    p2 <- build_schedule(p1, ctrl = ctrl)
    scm <- run_scm(p2, ctrl = ctrl)
    r_rr <- scm$net_reproduction_ratios

    scm$run_mutant(p2)
    m_rr <- scm$net_reproduction_ratios

    dplyr::tibble(birth_rate = x, resident_f = log(r_rr), mutant_f = log(m_rr))
  }

  p0$max_patch_lifetime <- 105.32 # default
  pr1 <- expand_parameters(traits, p0, birth_rate_list = 1)

  expected_eq <- 17.31739
  birth_rates <- c(0, 10, expected_eq * c(0.995, 1, 1.005), 25)

  outputs <- purrr::map_df(birth_rates, ~ f_test(pr1, .x, traits))

  expect_equal(birth_rates, outputs$birth_rate, tol = tol)
  expect_equal(outputs$resident_f, outputs$mutant_f, tol = tol)

  # Test 2
  p0$max_patch_lifetime <- 50
  expected_eq <- 1.662159
  birth_rates <- c(0, 1, expected_eq * c(0.995, 1, 1.005), 1.1)

  pr1 <- expand_parameters(traits, p0, birth_rate_list = 1)

  outputs <- purrr::map_df(birth_rates, ~ f_test(pr1, .x, traits))

  expect_equal(birth_rates, outputs$birth_rate, tol = tol)
  expect_equal(outputs$resident_f, outputs$mutant_f, tol = tol)
  
})
