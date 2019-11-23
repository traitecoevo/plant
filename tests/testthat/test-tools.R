context("Tools")
strategy_types <- get_list_of_strategy_types()

test_that("fixed_environment", {
  env <- fixed_environment(0.5)
  expect_is(env, "Environment")
    expect_equal(env$light_environment$xy, cbind(c(0, 75, 150), 0.5))
  expect_equal(env$canopy_openness(40), 0.5)
})

test_that("lcp_whole_plant", {
  for (x in names(strategy_types)) {
    ## R implementation:
    lcp_whole_plant_R <- function(plant, ...) {
      target <- function(canopy_openness) {
        env <- fixed_environment(canopy_openness)
        plant$compute_rates(env)
        plant$aux("net_mass_production_dt")
      }

      f1 <- target(1)
      if (f1 < 0.0) {
        NA_real_
      } else {
        uniroot(target, c(0, 1), f.upper=f1, ...)$root
      }
    }

    p <- Plant(x)(strategy_types[[x]]())
    expect_equal(lcp_whole_plant(p), lcp_whole_plant_R(p), tolerance=1e-5)
  }
})
