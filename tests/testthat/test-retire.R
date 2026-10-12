## Retiring cohorts whose density has collapsed (#519). A long run's cost is
## otherwise proportional to every cohort ever introduced, because a cohort with
## no density still pays a rate evaluation per stage; the pilot in
## notes/plan-519-continuous-patch.md measured one TF24 stand at 2.7 ms per
## introduced cohort per simulated year with 20 of 1000 cohorts carrying any.

ff16 <- function(lifetime = 60) {
  p <- scm_base_parameters("FF16")
  p$max_patch_lifetime <- lifetime
  p$node_schedule_times <- list()
  add_strategies(p, trait_matrix(0.0825, "lma"))
}

run_collecting <- function(p, ctrl) {
  scm <- SCM("FF16", "FF16_Env")(p, Environment("FF16"), empty_events(), ctrl)
  scm$collect <- TRUE
  scm$run()
  scm
}

test_that("retirement is off by default and needs the birth-date coordinate", {
  expect_equal(Control()$retire_survival, 0)
  p <- ff16(5)
  expect_error(run_scm(p, ctrl = Control(retire_survival = 1e-8)),
               "birth-date")
  expect_no_error(run_scm(p, ctrl = Control(retire_survival = 1e-8,
                                            node_density_in_birth_date = TRUE)))
})

test_that("cohorts below the floor are retired, frozen, and cost the stand less than the floor", {
  p <- ff16(60)
  ctrl0 <- Control(node_density_in_birth_date = TRUE)
  ctrl1 <- Control(node_density_in_birth_date = TRUE, retire_survival = 1e-6)
  scm0 <- run_collecting(p, ctrl0)
  scm1 <- run_collecting(p, ctrl1)
  sp0 <- scm0$patch$species[[1]]
  sp1 <- scm1$patch$species[[1]]

  ## Some cohorts went below the floor and were retired; none were without it.
  retired <- sp1$retired
  expect_gt(sum(retired), 0)
  expect_false(any(sp0$retired))
  ## A retired cohort carries exactly zero density, in the state the solver
  ## integrates, and its individual froze: its height stopped where it was.
  expect_true(all(exp(sp1$log_densities_state[retired]) == 0))
  ## Per retired cohort, over the snapshots in which it exists: once retired,
  ## its height never moves again.
  for (j in which(retired)) {
    exists <- vapply(scm1$history, function(hh) hh$species[[1]]$size >= j, logical(1))
    h <- vapply(scm1$history[exists], function(hh) hh$species[[1]]$heights[[j]], numeric(1))
    r <- vapply(scm1$history[exists], function(hh) hh$species[[1]]$retired[[j]], logical(1))
    expect_true(any(r))
    after <- which(r)
    expect_equal(h[after], rep(h[after[[1]]], length(after)))
  }
  ## The living cohorts are not frozen.
  expect_gt(max(abs(sp1$heights[!retired] - sp0$heights[!retired])), 0)
  ## What the stand loses is bounded by the floor: a retired cohort's remaining
  ## output is at most its remaining per-capita output times the floor over its
  ## density at birth, and the floor is six orders below the birth density here.
  expect_equal(scm1$net_reproduction_ratios, scm0$net_reproduction_ratios,
               tolerance = 1e-6)
  expect_equal(sum(exp(sp1$log_densities_state)), sum(exp(sp0$log_densities_state)),
               tolerance = 1e-6)
})

test_that("a retired cohort stays retired across an export and a resume", {
  p <- ff16(60)
  ctrl <- Control(node_density_in_birth_date = TRUE, retire_survival = 1e-6)
  scm <- run_collecting(p, ctrl)
  retired_end <- scm$patch$species[[1]]$retired
  ## A snapshot after the first retirement.
  first <- which(vapply(scm$history, function(h) any(h$species[[1]]$retired), logical(1)))[[1]]
  state <- export_patch_state(scm, step = first)
  p2 <- set_initial_state(p, state)
  scm2 <- SCM("FF16", "FF16_Env")(p2, Environment("FF16"), empty_events(), ctrl)
  scm2$collect <- TRUE
  scm2$run()
  ## The seeded patch recognises its retired cohorts from the state alone.
  expect_identical(scm2$history[[1]]$species[[1]]$retired,
                   scm$history[[first]]$species[[1]]$retired)
  ## And the resumed run ends with the same cohorts retired as the straight one.
  expect_identical(scm2$patch$species[[1]]$retired, retired_end)
  expect_equal(scm2$net_reproduction_ratios, scm$net_reproduction_ratios,
               tolerance = 1e-3)
})
