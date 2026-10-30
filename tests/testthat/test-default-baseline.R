## The continuous patch's features (#519) are opt-in, so a run that uses none of
## them must not move. Two checks: in-build invariances that hold on every
## platform, and a baseline recorded from the build before the features
## (scripts/capture-default-baseline.R), compared where it was recorded.

test_that("default runs reproduce the baseline recorded before #661", {
  skip_on_cran()
  ## Recorded on macOS arm64. Elsewhere the same code rounds differently and an
  ## adaptive step can land elsewhere, so the comparison would test the platform.
  skip_if_not(Sys.info()[["sysname"]] == "Darwin" && R.version$arch == "aarch64",
              "baseline recorded on macOS arm64")
  baseline <- readRDS(test_path("test_data", "default_baseline.rds"))
  current <- default_baseline_cases()
  expect_identical(names(current), names(baseline))
  for (case in names(baseline)) {
    for (what in names(baseline[[case]])) {
      expect_equal(current[[case]][[what]], baseline[[case]][[what]],
                   tolerance = 1e-12, label = paste(case, what))
    }
  }
})

test_that("snapshots read the run and do not move it", {
  run <- function(interval) {
    p <- scm_base_parameters("FF16")
    p$patch_type <- "fixed"
    p$max_patch_lifetime <- 60
    p <- add_strategies(p, trait_matrix(0.0825, "lma"), birth_rate = 1)
    p$node_schedule_times <- list(node_schedule_times_periodic(0, 60, 1))
    s <- run_scm(p, ctrl = Control(node_density_coordinate = "birth_date"), run = FALSE)
    s$snapshot_interval <- interval
    s$run()
    s
  }
  a <- run(0)
  b <- run(5)
  expect_gt(length(b$snapshots), 10)
  expect_identical(b$patch$ode_state, a$patch$ode_state)
  expect_identical(b$ode_times, a$ode_times)
  expect_identical(b$net_reproduction_ratios, a$net_reproduction_ratios)
})

test_that("a fitness window weights the run and does not move it", {
  ## The window is a weighting of arrivals: the stand it weights is the stand
  ## an unweighted run makes, step for step.
  p <- scm_base_parameters("FF16")
  p$patch_type <- "fixed"
  p$max_patch_lifetime <- 60
  p <- add_strategies(p, trait_matrix(0.0825, "lma"), birth_rate = 1)
  p$node_schedule_times <- list(node_schedule_times_periodic(0, 60, 1))
  ctrl <- Control(node_density_coordinate = "birth_date")
  a <- run_scm(p, ctrl = ctrl)
  b <- run_scm(fitness_window(p, 10, 30), ctrl = ctrl)
  expect_identical(b$patch$ode_state, a$patch$ode_state)
  expect_identical(b$ode_times, a$ode_times)
  expect_false(identical(b$net_reproduction_ratios, a$net_reproduction_ratios))
})
