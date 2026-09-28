# stand_gradient(hyperpar = ): the total derivative with respect to a trait as a
# trait matrix gives it, through the parameters the hyperparameter function
# derives from it. Without `hyperpar` a column is a partial with every derived
# parameter held fixed, and for lma that was 2.4 to 4 times the derivative a
# trait fit wants.

hyperpar_run <- function(lifetime = 3) {
  p <- scm_base_parameters("TF24")
  p$max_patch_lifetime <- lifetime
  p <- add_strategies(p, trait_matrix(c(0.0825, 5.13), c("lma", "hmat")),
                      hyperpar = TF24_hyperpar, birth_rate = list(1.10))
  ctrl <- Control(node_density_in_birth_date = TRUE)
  list(p = p, ctrl = ctrl,
       scm = run_scm(p, Environment("TF24"), ctrl, refine_schedule = TRUE,
                     record_trajectory = TRUE))
}

test_that("the chain rule is the partials weighted by the hyperpar's own derivative", {
  run <- hyperpar_run()
  tot <- stand_gradient(run$scm, traits = c("1.lma", "1.hmat"),
                        hyperpar = TF24_hyperpar)
  par <- stand_gradient(run$scm)

  # lma's derived parameters, and one checked against TF24_hyperpar's closed
  # form: k_l = B_kl1 * (lma / lma_0)^(-B_kl2), so dk_l/dlma = -B_kl2 * k_l / lma.
  J <- tot$jacobian[["1.lma"]]
  expect_true(all(c("k_l", "r_l") %in% names(J)))
  lma <- 0.0825
  k_l <- run$scm$parameters$strategies[[1]]$pars$k_l
  expect_equal(J[["k_l"]], -1.71 * k_l / lma, tolerance = 1e-8)

  # The total, rebuilt from the partials: own column plus each derived column
  # times its derivative. Parameters no equation reads carry no column.
  cols <- paste0("1.", names(J))
  keep <- cols %in% colnames(par$gradient)
  expect_equal(tot$gradient[, "1.lma"],
               par$gradient[, "1.lma"] +
                 as.vector(par$gradient[, cols[keep], drop = FALSE] %*% J[keep]))

  # Nothing derives from hmat, so its total is its partial.
  expect_length(tot$jacobian[["1.hmat"]], 0L)
  expect_identical(tot$gradient[, "1.hmat"], par$gradient[, "1.hmat"])

  # And the point of it: for lma the partial is not the total.
  expect_true(all(abs(par$gradient[, "1.lma"] / tot$gradient[, "1.lma"]) > 2))
})

test_that("the total agrees with a difference of the trait through add_strategies", {
  skip_on_cran()
  run <- hyperpar_run()
  tot <- stand_gradient(run$scm, traits = "1.lma", hyperpar = TF24_hyperpar)
  sched <- run$scm$node_schedule$times(1)
  census_at <- function(lma) {
    p <- scm_base_parameters("TF24")
    p$max_patch_lifetime <- 3
    p <- add_strategies(p, trait_matrix(c(lma, 5.13), c("lma", "hmat")),
                        hyperpar = TF24_hyperpar, birth_rate = list(1.10))
    p$node_schedule_times <- list(sched)
    stand_census(run_scm(p, Environment("TF24"), run$ctrl,
                         refine_schedule = FALSE))
  }
  h <- 1e-3 * 0.0825
  fd <- (census_at(0.0825 + h) - census_at(0.0825 - h)) / (2 * h)
  # Measured 3e-5 to 8e-5 across two steps, which is the difference's own noise;
  # the partial it replaces was 2.4 to 3.4 times off.
  expect_equal(tot$gradient[names(fd), "1.lma"], fd, tolerance = 1e-3)
})

test_that("a trait it cannot differentiate through is refused", {
  run <- hyperpar_run(lifetime = 1)
  expect_error(stand_gradient(run$scm, hyperpar = TF24_hyperpar),
               "Name the traits")
  expect_error(stand_gradient(run$scm, traits = "lma", hyperpar = TF24_hyperpar),
               "names its species")
  expect_error(stand_gradient(run$scm, traits = "1.lma", hyperpar = "TF24"),
               "must be a hyperparameter function")
})
