if (interactive()) {
  devtools::load_all("../../")
  library(testthat)
  source("helper-tree2.R")
}

context("PlantRunner")

test_that("PlantRunner", {
  p <- Plant(Strategy())
  env <- test_environment(10)

  pr <- PlantRunner(p, env)
  expect_that(pr, is_a("PlantRunner"))
  expect_that(pr$plant, is_a("Plant"))
  expect_that(pr$plant$vars_size, is_identical_to(p$vars_size))

  ## This going to work with a *copy* of pr; so that won't propagate
  ## back.
  runner <- OdeRunner("PlantRunner")(pr)
  expect_that(runner, is_a("OdeRunner"))
  expect_that(runner, is_a("OdeRunner<PlantRunner>"))
  expect_that(runner$time, equals(0.0))


  continue_if <- function(obj) {
    obj$state[[1]] < 15
  }
  observer <- function(obj) {
    c(obj$time, obj$state)
  }
  pr <- PlantRunner(Plant(Strategy()), env)
  runner <- OdeRunner("PlantRunner")(pr)
  ret <- list(observer(runner))
  while (continue_if(runner)) {
    message(runner$time)
    runner$step()
    ret <- c(ret, list(observer(runner)))
  }

  m <- do.call("rbind", ret)
  colnames(m) <- c("time", runner$object$plant$ode_names)

  if (interactive()) {
    plot(height ~ time, as.data.frame(m))
  }
})
