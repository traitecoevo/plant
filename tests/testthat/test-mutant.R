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
  birth_rate <- 100
  
  pn <- expand_parameters(trait_matrix(rep(lma, n), "lma"), p,
                          mutant = F,
                          birth_rate_list = rep(birth_rate, n))
  
  types <- extract_RcppR6_template_types(pn, "Parameters")
  scm <- do.call('SCM', types)(pn, e, ctrl)
  scm$run()
  
  # check fitness
  resident_rr <- scm$net_reproduction_ratios
  expect_equal(resident_rr, c(4.290414e-23, 4.290414e-23))
  
  
  # identical mutants - same environment
  scm$run_mutant(pn)
  expect_setequal(resident_rr, scm$net_reproduction_ratios)
  
  # one mutant
  p1 <- expand_parameters(trait_matrix(lma, "lma"), p,
                          mutant = F,
                          birth_rate_list = birth_rate)
  
  scm$run_mutant(p1)
  expect_equal(resident_rr[1], scm$net_reproduction_ratios)
  
  # n + 1 mutants
  pn_p1 <- expand_parameters(trait_matrix(c(rep(lma, n + 1)), "lma"), p,
                             mutant = F,
                             birth_rate_list = rep(birth_rate, n + 1))
  
  scm$run_mutant(pn_p1)
  expect_equal(rep(resident_rr[1], n + 1), scm$net_reproduction_ratios)
  
  # non-identical mutant
  new_lma = 0.20
    
  p1 <- expand_parameters(trait_matrix(new_lma, "lma"), p,
                          mutant = F,
                          birth_rate_list = birth_rate)
  
  scm$run_mutant(p1)
  expect_false(resident_rr[1] == scm$net_reproduction_ratios)
})
