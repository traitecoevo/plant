# A trait's column is the total derivative as the trait moves through the
# hyperparameter function that built the strategy; hyperpar = NULL gives the
# partials, with every derived parameter fixed. For lma the partial is 2.4
# to 4 times the total.

hyperpar_parameters <- function(lifetime = 3, lma = 0.0825) {
  p <- testing_parameters("TF24")
  p$max_patch_lifetime <- lifetime
  add_strategies(p, trait_matrix(c(lma, 5.13), c("lma", "hmat")),
                 hyperpar = default_hyperpar("TF24"), birth_rate = list(1.10))
}

hyperpar_run <- function(p = hyperpar_parameters()) {
  ctrl <- Control(node_density_coordinate = "birth_date")
  list(p = p, ctrl = ctrl,
       solver = run(p, Environment("TF24"), ctrl, refine_schedule = TRUE,
                    record = "states")$solver)
}

# The default run, made once: two tests read it, and neither changes it (a
# sweep reads a finished run's record and leaves it as it was).
default_hyperpar_run <- local({
  run <- NULL
  function() {
    if (is.null(run)) run <<- hyperpar_run()
    run
  }
})

test_that("a trait's column is its partial plus each derived parameter's, chained", {
  run <- default_hyperpar_run()
  tot <- stand_gradient(run$solver, traits = c("1.lma", "1.hmat"))
  par <- stand_gradient(run$solver, hyperpar = NULL)

  # the TF24 trade-off function derives k_l and r_l from lma. Their derivatives,
  # from the function itself.
  s <- run$solver$parameters$strategies[[1]]
  derived <- function(lma) {
    default_hyperpar("TF24")(trait_matrix(lma, "lma"), s, filter = FALSE)[1, c("k_l", "r_l")]
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
  run <- default_hyperpar_run()
  tot <- stand_gradient(run$solver, traits = "1.lma")
  # The trait moves through add_strategies, on the run's own steps.
  pinned <- pin_steps(run$solver)
  census_at <- function(lma) {
    p <- pinned
    p$strategies <- hyperpar_parameters(lma = lma)$strategies
    stand_census(run(p, Environment("TF24"), run$ctrl)$solver)
  }
  h <- 1e-4 * 0.0825
  fd <- (census_at(0.0825 + h) - census_at(0.0825 - h)) / (2 * h)
  # Measured 3e-7 to 4e-7; the partial is 2.4 to 3.4 times off.
  for (m in names(fd)) {
    expect_equal(tot$gradient[m, "1.lma"], fd[[m]], tolerance = 1e-5, label = m)
  }
})

test_that("a derived parameter is not a column of its own under a trade-off", {
  run <- hyperpar_run(hyperpar_parameters(lifetime = 1))
  tot <- stand_gradient(run$solver)
  par <- stand_gradient(run$solver, hyperpar = NULL)
  expect_false(any(c("1.k_l", "1.r_l", "1.stem_P50") %in% colnames(tot$gradient)))
  expect_true(all(c("1.lma", "1.K_s", "1.hmat") %in% colnames(tot$gradient)))
  expect_true(all(c("1.k_l", "1.r_l") %in% colnames(par$gradient)))
  expect_error(stand_gradient(run$solver, traits = "1.k_l"),
               "1.k_l derived from the species' traits")
})

test_that("a trade-off deriving a parameter that is read and has no column is refused", {
  run <- hyperpar_run(hyperpar_parameters(lifetime = 1))
  # The same trade-off, also deriving one more parameter at the value the
  # strategy holds, so the strategy is on it.
  deriving <- function(name) {
    hp <- default_hyperpar("TF24")
    function(m, s, filter = TRUE) {
      out <- hp(m, s, filter = filter)
      out <- cbind(out, rep(s$pars[[name]], nrow(out)))
      colnames(out)[ncol(out)] <- name
      out
    }
  }
  # g_min reaches the leaf and has no column, so its channel would be lost.
  expect_error(stand_gradient(run$solver, traits = "1.lma",
                              hyperpar = deriving("g_min")),
               "derives 1.g_min .*chain rule cannot be completed")
  # root_psi_crit reaches no equation, so the chain loses nothing by it.
  # (psi_crit would test nothing: the default trade-off derives it already.)
  expect_no_error(stand_gradient(run$solver, traits = "1.lma",
                                 hyperpar = deriving("root_psi_crit")))
})

test_that("a strategy that is not on the trade-off is refused, not chained", {
  p <- hyperpar_parameters(lifetime = 1)
  p <- add_strategies(p, trait_matrix(c(0.0825, 5.13), c("lma", "hmat")),
                      hyperpar = default_hyperpar("TF24"), birth_rate = list(1.10))
  p$strategies[[2]]$pars$lma <- 0.09
  run <- hyperpar_run(p)
  expect_error(stand_gradient(run$solver, traits = "2.lma"),
               "Species 2's k_l, r_l are not what `hyperpar` gives")
  expect_error(stand_gradient(run$solver), "Species 2's")
  # Only the species asked about is checked.
  expect_no_error(stand_gradient(run$solver, traits = "1.lma"))
  expect_no_error(stand_gradient(run$solver, traits = "2.lma", hyperpar = NULL))
  expect_error(stand_gradient(run$solver, traits = "1.lma", hyperpar = "TF24"),
               "must be a trade-off function")
})

test_that("a default strategy, or one built at an input it does not keep, is refused", {
  # TF24's defaults are not the trade-off's at its default traits.
  expect_error(trade_off(testing_parameters("TF24")$strategy_default, 1,
                         default_hyperpar("TF24")), "built by hand")
  p <- add_strategies(testing_parameters("FF16"),
                      trait_matrix(c(0.1, 3e-3), c("lma", "narea")),
                      hyperpar = default_hyperpar("FF16"), birth_rate = list(1))
  expect_error(trade_off(p$strategies[[1]], 1, default_hyperpar("FF16")),
               "built at a narea the strategy does not keep")
  expect_no_error(trade_off(p$strategies[[1]], 1, default_hyperpar("FF16", narea = 3e-3)))
})
