# A trait's column is the total derivative as the trait moves through the
# hyperparameter function that built the strategy; hyperpar = NULL gives the
# partials, with every derived parameter fixed. For lma the partial was 2.4
# to 4 times the total.

hyperpar_parameters <- function(lifetime = 3) {
  p <- scm_base_parameters("TF24")
  p$max_patch_lifetime <- lifetime
  add_strategies(p, trait_matrix(c(0.0825, 5.13), c("lma", "hmat")),
                 hyperpar = TF24_hyperpar, birth_rate = list(1.10))
}

hyperpar_run <- function(p = hyperpar_parameters()) {
  ctrl <- Control(node_density_coordinate = "birth_date")
  list(p = p, ctrl = ctrl,
       scm = run_scm(p, Environment("TF24"), ctrl, refine_schedule = TRUE,
                     keep_states = TRUE))
}

test_that("a trait's column is its partial plus each derived parameter's, chained", {
  run <- hyperpar_run()
  tot <- stand_gradient(run$scm, traits = c("1.lma", "1.hmat"))
  par <- stand_gradient(run$scm, hyperpar = NULL)

  # TF24_hyperpar derives k_l, r_l and nmass_l from lma; nmass_l reaches no
  # equation, so it has no column. Their derivatives, from the function itself.
  s <- run$scm$parameters$strategies[[1]]
  derived <- function(lma) {
    TF24_hyperpar(trait_matrix(lma, "lma"), s, filter = FALSE)[1, c("k_l", "r_l")]
  }
  h <- 1e-6 * 0.0825
  d <- (derived(0.0825 + h) - derived(0.0825 - h)) / (2 * h)
  # k_l = B_kl1 * (lma / lma_0)^(-B_kl2), so dk_l/dlma = -B_kl2 * k_l / lma.
  expect_equal(d[["k_l"]], -1.71 * s$pars$k_l / 0.0825, tolerance = 1e-8)
  expect_equal(tot$gradient[, "1.lma"],
               par$gradient[, "1.lma"] + par$gradient[, "1.k_l"] * d[["k_l"]] +
                 par$gradient[, "1.r_l"] * d[["r_l"]], tolerance = 1e-8)

  # Nothing derives from hmat, so its total is its partial.
  expect_identical(tot$gradient[, "1.hmat"], par$gradient[, "1.hmat"])

  # And the point of it: for lma the partial is not the total.
  expect_true(all(abs(par$gradient[, "1.lma"] / tot$gradient[, "1.lma"]) > 2))
})

test_that("the total agrees with a difference of the trait through add_strategies", {
  skip_on_cran()
  run <- hyperpar_run()
  tot <- stand_gradient(run$scm, traits = "1.lma")
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
  # the partial was 2.4 to 3.4 times off.
  expect_equal(tot$gradient[names(fd), "1.lma"], fd, tolerance = 1e-3)
})

test_that("a derived parameter is not a column of its own under a trade-off", {
  run <- hyperpar_run(hyperpar_parameters(lifetime = 1))
  tot <- stand_gradient(run$scm)
  par <- stand_gradient(run$scm, hyperpar = NULL)
  expect_false(any(c("1.k_l", "1.r_l", "1.stem_P50") %in% colnames(tot$gradient)))
  expect_true(all(c("1.lma", "1.K_s", "1.hmat") %in% colnames(tot$gradient)))
  expect_true(all(c("1.k_l", "1.r_l") %in% colnames(par$gradient)))
  expect_error(stand_gradient(run$scm, traits = "1.k_l"),
               "1.k_l derived from the species' traits")
})

test_that("a strategy that is not on the trade-off is refused, not chained", {
  p <- hyperpar_parameters(lifetime = 1)
  p <- add_strategies(p, trait_matrix(c(0.0825, 5.13), c("lma", "hmat")),
                      hyperpar = TF24_hyperpar, birth_rate = list(1.10))
  p$strategies[[2]]$pars$lma <- 0.09
  run <- hyperpar_run(p)
  expect_error(stand_gradient(run$scm, traits = "2.lma"),
               "Species 2's k_l, r_l, nmass_l are not what `hyperpar` gives")
  expect_error(stand_gradient(run$scm), "Species 2's")
  # Only the species asked about is checked.
  expect_no_error(stand_gradient(run$scm, traits = "1.lma"))
  expect_no_error(stand_gradient(run$scm, traits = "2.lma", hyperpar = NULL))
  expect_error(stand_gradient(run$scm, traits = "1.lma", hyperpar = "TF24"),
               "must be a hyperparameter function")
})

test_that("a default strategy, or one built at an input it does not keep, is refused", {
  # TF24's defaults are not the trade-off's at its default traits.
  expect_error(trade_off(scm_base_parameters("TF24")$strategy_default, 1,
                         TF24_hyperpar), "built by hand")
  p <- add_strategies(scm_base_parameters("FF16"),
                      trait_matrix(c(0.1, 3e-3), c("lma", "narea")),
                      hyperpar = FF16_hyperpar, birth_rate = list(1))
  expect_error(trade_off(p$strategies[[1]], 1, FF16_hyperpar),
               "built at a narea the strategy does not keep")
  expect_no_error(trade_off(p$strategies[[1]], 1, make_FF16_hyperpar(narea = 3e-3)))
})
