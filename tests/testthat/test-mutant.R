test_that("mutant runs work", {
  # basic setup 
  p <- scm_base_parameters("FF16")
  p$max_patch_lifetime <- 50
  
  e <- make_environment("FF16")
  ctrl <- scm_base_control()
  ctrl$save_RK45_cache = TRUE
    
  tol <- 1e-4
  
  # 3 residents

  lma <- c(0.05, 0.1, 0.2)
  birth_rate <- 1

  # 1 resident strategies
  p1 <- expand_parameters(trait_matrix(lma[2], "lma"), p,
    birth_rate_list = rep(birth_rate, 1)
  )

  p1m1 <- expand_parameters(trait_matrix(lma[3], "lma"), p1,
    birth_rate_list = rep(birth_rate, 1)
  )

  p1m3 <- expand_parameters(trait_matrix(lma, "lma"), p1,
    birth_rate_list = rep(birth_rate, 3)
  )

  p1m10 <- expand_parameters(trait_matrix(seq(lma[1], lma[3], length.out=10), "lma"), p1,
    birth_rate_list = rep(birth_rate, 10)
  )

  types <- extract_RcppR6_template_types(p1, "Parameters")
  scm <- do.call("SCM", types)(p1, e, ctrl)

  # test error handling
  expect_error(scm$run_mutant(p), "Run a resident first to generate a competitve landscape") 

  scm$run()
  p1_rr <- scm$net_reproduction_ratios
  expected <- 2.77322
  expect_equal(p1_rr, expected, tol = tol)

  scm$run_mutant(p1m1)
  p1m1_rr <- scm$net_reproduction_ratios
  expected <- c(2.77322, 3.707605)
  expect_equal(p1m1_rr, expected, tol = tol)
  expect_equal(p1m1_rr[1], p1_rr, tol = tol)

  scm$run_mutant(p1m3)
  p1m3_rr <- scm$net_reproduction_ratios
  expected <- c(2.77322, 3.7429e-10, 2.77322, 3.70753)
  expect_equal(p1m3_rr, expected, tol = tol)
  expect_equal(p1m3_rr[1], p1_rr, tol = tol)

  scm$run_mutant(p1m10)
  p1m10_rr <- scm$net_reproduction_ratios
  expected <- c(2.773222, 3.742935e-10, 9.308944e-07, 0.1363641, 2.773222, 3.890554, 1.524582, 1.160212, 1.871261, 2.765328, 3.707372)
  expect_equal(p1m10_rr, expected, tol = tol)
  expect_equal(p1m10_rr[1], p1_rr, tol = tol)

  # 3 resident strategies
  p3 <- expand_parameters(trait_matrix(lma, "lma"), p,
    birth_rate_list = rep(birth_rate, 3)
  )
  
  p3m1 <- expand_parameters(trait_matrix(lma[3], "lma"), p3,
    birth_rate_list = rep(birth_rate, 1)
  )

  p3m3 <- expand_parameters(trait_matrix(lma, "lma"), p3,
    birth_rate_list = rep(birth_rate, 3)
  )

  p3m10 <- expand_parameters(trait_matrix(seq(lma[1], lma[3], length.out = 10), "lma"), p3,
    birth_rate_list = rep(birth_rate, 10)
  )

  types <- extract_RcppR6_template_types(p3, "Parameters")
  scm <- do.call("SCM", types)(p3, e, ctrl)
  scm$run()
  p3_rr <- scm$net_reproduction_ratios
  expected <- c(4.265e-10, 2.831741, 0.09125339)
  expect_equal(p3_rr, expected, tol = tol)


  scm$run_mutant(p3m1)
  p3m1_rr <- scm$net_reproduction_ratios
  expected <- c(4.265e-10, 2.831741, 0.09125339, 0.09125339)
  expect_equal(p3m1_rr, expected, tol = tol)
  expect_equal(p3m1_rr[1:3], p3_rr, tol = tol)

  scm$run_mutant(p3m3)
  p3m3_rr <- scm$net_reproduction_ratios
  expected <- c(4.265e-10, 2.831741, 0.09125339, 4.265e-10, 2.831741, 0.09125339)
  expect_equal(p3m3_rr, expected, tol = tol)
  expect_equal(p3m3_rr[1:3], p3_rr, tol = tol)

  scm$run_mutant(p3m10)
  p3m10_rr <- scm$net_reproduction_ratios
  expected <- c(4.265011e-10, 2.831741, 0.09125377, 4.265011e-10, 5.587752e-06, 0.266188, 2.831741, 2.690585, 0.3796333, 0.07098642, 0.07226859, 0.08342181, 0.09125377)
  expect_equal(p3m10_rr, expected, tol = tol)
  expect_equal(p3m3_rr[1:3], p3_rr, tol = tol)

  # Test with different resident densities
  
  
  # 1 resident strategies
  ctrl <- scm_base_control()
  ctrl$equilibrium_eps <- 1e-5
  ctrl$save_RK45_cache = TRUE


  f <- function(p, x, traits) {
    p2 <- p
    p2$strategies[[1]]$birth_rate_y <- x

    p3 <- build_schedule(p2, ctrl = ctrl)
    scm <- run_scm(p3, ctrl = ctrl)
    r_rr <- scm$net_reproduction_ratios

    scm$run_mutant(p3)
    m_rr <- scm$net_reproduction_ratios

    fl1 <- fitness_landscape(traits, p2)
    fl2 <- fitness_landscape(traits, p3)
    dplyr::tibble(birth_rate = x, resident_f = log(r_rr), mutant_f = log(m_rr), landscape_f1 = fl1, landscape_f2 = fl2)
  } 

  lma_attr <- 0.0825
  p <- scm_base_parameters("FF16")
  traits <- trait_matrix(lma_attr, c("lma"))

  p$max_patch_lifetime <- 105.32 # default

  p1 <- expand_parameters(traits, p, birth_rate_list = 1)

  expected_eq <- 17.31739
  birth_rates <- c(seq(0, 25, length.out = 6), expected_eq * c(0.995, 1, 1.005)) %>% sort()

  p_eq <- equilibrium_birth_rate(p1, ctrl=ctrl)

  expect_equal(
    p_eq$strategies[[1]]$birth_rate_y,
    expected_eq, tol=tol)

  outputs <- purrr::map_df(birth_rates, ~ f(p_eq, .x, traits))

  outputs
  
  expect_equal(birth_rates, outputs$birth_rate, tol=tol)
  expect_equal(outputs$resident_f, outputs$mutant_f, tol = tol * 10)
  expect_equal(outputs$resident_f, outputs$landscape_f1, tol = tol * 10)
  expect_equal(outputs$resident_f, outputs$landscape_f2, tol = tol * 10)

  # Test 2

  p$max_patch_lifetime <- 50
  expected_eq <- 1.662159
  birth_rates <- c(seq(0, 2, length.out = 6), expected_eq * c(0.995, 1, 1.005)) %>% sort()

  p1 <- expand_parameters(traits, p, birth_rate_list = 1)

  p_eq <- equilibrium_birth_rate(p1, ctrl=ctrl)

  expect_equal(
    p_eq$strategies[[1]]$birth_rate_y,
    expected_eq, tol=tol)

  outputs <- purrr::map_df(birth_rates, ~ f(p_eq, .x, traits))

  expect_equal(birth_rates, outputs$birth_rate, tol=tol)
  expect_equal(outputs$resident_f, outputs$mutant_f, tol = tol * 10)
  expect_equal(outputs$resident_f, outputs$landscape_f1, tol = tol * 10)
  expect_equal(outputs$resident_f, outputs$landscape_f2, tol = tol * 10)


  # Test with equilbirum seed rain

  #XXXX

  # Test with build_schedule

  #XXXX


})
