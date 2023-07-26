test_that("mutant runs work", {
  # basic setup 
  p <- scm_base_parameters("FF16")
  p$max_patch_lifetime <- 1
  
  e <- make_environment("FF16")
  ctrl <- scm_base_control()
  ctrl$save_RK45_cache = TRUE
    
  
  # n identical residents
  lma <- 0.19
  n <- 2
  birth_rate <- 1
  
  pn <- expand_parameters(trait_matrix(rep(lma, n), "lma"), p,
                          mutant = F,
                          birth_rate_list = rep(birth_rate, n))
  
  types <- extract_RcppR6_template_types(pn, "Parameters")
  scm <- do.call('SCM', types)(pn, e, ctrl)
  scm$run()
  
  # check fitness
  resident_rr <- scm$net_reproduction_ratios
  expect_equal(resident_rr, c(4.520153e-23, 4.520153e-23), tolerance = 1e-24)
  
  
  # identical mutants - same environment
  scm$run_mutant(pn)
  expect_equal(resident_rr, scm$net_reproduction_ratios, tolerance = 1e-24)
  
  # one mutant
  p1 <- expand_parameters(trait_matrix(lma, "lma"), p,
                          mutant = F,
                          birth_rate_list = birth_rate)
  
  scm$run_mutant(p1)
  expect_equal(resident_rr[1], scm$net_reproduction_ratios, tolerance = 1e-24)
  
  # n + 1 mutants
  pn_p1 <- expand_parameters(trait_matrix(c(rep(lma, n + 1)), "lma"), p,
                             mutant = F,
                             birth_rate_list = rep(birth_rate, n + 1))
  
  scm$run_mutant(pn_p1)
  expect_equal(rep(resident_rr[1], n + 1), scm$net_reproduction_ratios, tolerance = 1e-24)
  
  # non-identical mutant - non-invasible
  new_lma = 0.20
    
  p1 <- expand_parameters(trait_matrix(new_lma, "lma"), p,
                          mutant = F,
                          birth_rate_list = birth_rate)
  
  scm$run_mutant(p1)
  expect_gt(resident_rr[1] - scm$net_reproduction_ratios, 0.0)
  
  # non-identical mutant - viable invader
  new_lma = 0.18
  
  p1 <- expand_parameters(trait_matrix(new_lma, "lma"), p,
                          mutant = F,
                          birth_rate_list = birth_rate)
  
  scm$run_mutant(p1)
  expect_lt(resident_rr[1] - scm$net_reproduction_ratios, 0)
})

test_that("many mutants for speed", {
  # basic setup
  p <- scm_base_parameters("FF16")
  p$max_patch_lifetime <- 100

  e <- make_environment("FF16")
  ctrl <- scm_base_control()
  ctrl$save_RK45_cache <- TRUE
  ctrl_no_cache <- scm_base_control()


  # n identical residents
  lma <- c(0.05, 0.1, 0.2)
  birth_rate <- 1

  p1 <- expand_parameters(trait_matrix(lma[2], "lma"), p,
    mutant = F,
    birth_rate_list = rep(birth_rate, 1)
  )

  p3 <- expand_parameters(trait_matrix(lma, "lma"), p,
    mutant = F,
    birth_rate_list = rep(birth_rate, 3)
  )

  p21 <- expand_parameters(trait_matrix(rep(lma, 7), "lma"), p,
    mutant = F,
    birth_rate_list = rep(birth_rate, 21)
  )

  types <- extract_RcppR6_template_types(p1, "Parameters")
  scm <- do.call("SCM", types)(p1, e, ctrl)
  cat("Residents 1: ", system.time(scm$run()), "\n")
  resident_rr1 <- scm$net_reproduction_ratios

  cat("mutants 1: ", system.time(scm$run_mutant(p1)), "\n")
  expect_equal(resident_rr1, scm$net_reproduction_ratios, tolerance = 1e-3)
  cat("mutants 3: ", system.time(scm$run_mutant(p3)), "\n")
  cat("mutants 21: ", system.time(scm$run_mutant(p21)), "\n")
 

  types <- extract_RcppR6_template_types(p3, "Parameters")
  scm <- do.call("SCM", types)(p3, e, ctrl)
  cat("Residents 3: ", system.time(scm$run()), "\n")
  resident_rr3 <- scm$net_reproduction_ratios

  cat("mutants 1: ", system.time(scm$run_mutant(p1)), "\n")
  cat("mutants 3: ", system.time(scm$run_mutant(p3)), "\n")
  expect_equal(resident_rr3, scm$net_reproduction_ratios, tolerance = 1e-3)
  cat("mutants 21: ", system.time(scm$run_mutant(p21)), "\n")


  types <- extract_RcppR6_template_types(p21, "Parameters")
  scm <- do.call("SCM", types)(p21, e, ctrl)
  cat("Residents 21: ", system.time(scm$run()), "\n")
  resident_rr21 <- scm$net_reproduction_ratios

  # one mutants - same environment
  cat("mutants 1: ", system.time(scm$run_mutant(p1)), "\n")
  cat("mutants 3: ", system.time(scm$run_mutant(p3)), "\n")
  cat("mutants 21: ", system.time(scm$run_mutant(p21)), "\n")
  expect_equal(resident_rr21, scm$net_reproduction_ratios, tolerance = 1e-3)

  scm <- do.call("SCM", types)(p1, e, ctrl_no_cache)
  cat("Residents 1: ", system.time(scm$run()), "\n")
  scm <- do.call("SCM", types)(p3, e, ctrl_no_cache)
  cat("Residents 3: ", system.time(scm$run()), "\n")
  scm <- do.call("SCM", types)(p21, e, ctrl_no_cache)
  cat("Residents 21: ", system.time(scm$run()), "\n")

})

