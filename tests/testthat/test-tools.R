context("Tools")
## TODO: Remove ["FFW16"] to test this with all types. But first
## requires issue #162 to be resolved
strategy_types <- get_list_of_strategy_types()

test_that("fixed_environment", {
  env <- fixed_environment(0.5)
  expect_that(env, is_a("Environment"))
    expect_that(env$light_environment$xy,
                equals(cbind(c(0, 75, 150), 0.5)))
  expect_that(env$canopy_openness(40), equals(0.5))
})

test_that("lcp_whole_plant", {
  for (x in names(strategy_types)) {
    ## R implementation:
    lcp_whole_plant_R <- function(plant, ...) {
      target <- function(canopy_openness) {
        env <- fixed_environment(canopy_openness)
        plant$compute_vars_phys(env)
        plant$internals[["net_mass_production_dt"]]
      }

      f1 <- target(1)
      if (f1 < 0.0) {
        NA_real_
      } else {
        uniroot(target, c(0, 1), f.upper=f1, ...)$root
      }
    }

    p <- PlantPlus(x)(strategy_types[[x]]())
    expect_that(lcp_whole_plant(p),
                equals(lcp_whole_plant_R(p), tolerance=1e-5))
  }
})
