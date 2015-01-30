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

  expect_that(oderunner_plant_size(runner),
              is_identical_to(p$vars_size))

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

test_that("grow_plant_to_size", {
  env <- test_environment(10)
  heights <- seq(1, 10)
  s <- Strategy()

  res <- grow_plant_bracket(Plant(s), heights, "height", env)

  expect_that(res$t0, is_identical_to(res$time[res$index]))
  expect_that(res$t1, is_identical_to(res$time[res$index + 1L]))
  expect_that(res$y0, is_identical_to(res$state[res$index,]))
  expect_that(res$y1, is_identical_to(res$state[res$index + 1L,]))
  ## We really do bracket the size:
  expect_that(all(res$y0[,"height"] < heights), is_true())
  expect_that(all(res$y1[,"height"] > heights), is_true())
  expect_that(res$runner, is_a("OdeRunner<PlantRunner>"))

  ## Then, do the search for a single case:
  i <- 3L
  tmp <- grow_plant_bisect(res$runner,
                           heights[[i]], "height",
                           res$t0[[i]], res$t1[[i]], res$y0[i,])

  ## The plant lies within the range expected:
  expect_that(tmp$time, is_within_interval(res$t0[[i]], res$t1[[i]]))
  expect_that(all(tmp$state > res$y0[i,]), is_true())
  expect_that(all(tmp$state < res$y1[i,]), is_true())

  ## Do all plants using the proper function:
  obj <- grow_plant_to_size(Plant(s), heights, "height", env)
  expect_that(obj$time, is_a("numeric"))
  expect_that(all(obj$time > res$t0), is_true())
  expect_that(all(obj$time < res$t1), is_true())

  expect_that(all(obj$state > res$y0), is_true())
  expect_that(all(obj$state < res$y1), is_true())

  expect_that(length(obj$plant), equals(length(heights)))
  expect_that(all(sapply(obj$plant, inherits, "Plant")), is_true())
  expect_that(sapply(obj$plant, function(p) p$height),
              equals(heights, tolerance=1e-6))
})

## TODO: another useful function could be to construct splines for
## arbitrary variables during a run.
