context("equilibrium birth rates")
test_that("demographic equilibrium", {
  params <- scm_base_parameters("FF16")
  ctrl <- scm_base_control()
  ctrl2 <- scm_base_control()
  
  patch <- expand_parameters(trait_matrix(0.0825, "lma"), params)
  expect_silent(
    patch_eq1 <- equilibrium_birth_rate(patch, ctrl=ctrl)
  )
  expect_true(attr(patch_eq1, "converged"))
  expect_equal(attr(patch_eq1, "offspring_production"), 17.31629, tolerance = 1e-4)
  
  ctrl2$equilibrium_solver_name <- "nleqslv"
  expect_silent(
    patch_eq2 <- equilibrium_birth_rate(patch, ctrl2)
  )
  expect_true(attr(patch_eq2, "converged"))
  expect_equal(attr(patch_eq2, "offspring_production"), 17.31729, tolerance = 1e-4)
  
  ctrl2$equilibrium_solver_name <- "hybrid"
  expect_silent(
    patch_eq3 <- equilibrium_birth_rate(patch, ctrl2)
  )
  expect_true(attr(patch_eq3, "converged"))
  expect_equal(attr(patch_eq3, "offspring_production"), 17.31727, tolerance = 1e-4)
  
  ctrl2$equilibrium_solver_name <- "nonsense"
  expect_error(
    patch_eq4 <- equilibrium_birth_rate(patch, ctrl2)
  )
  
  ## 2spp
  patch <- expand_parameters(trait_matrix(c(0.0825, 0.2), "lma"), params, birth_rate_list = c(1,1))
  expect_silent(
    patch_eq5 <- equilibrium_birth_rate(patch, ctrl)
  )
  expect_true(attr(patch_eq5, "converged"))
  expect_equal(attr(patch_eq5, "offspring_production"), c(12.41250, 13.9535), tolerance = 1e-4)
})