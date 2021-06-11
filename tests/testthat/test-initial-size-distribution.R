context("Initial size distribution")

# using SCM collect to generate starting conditions
p0 <- scm_base_parameters("FF16")
p0$patch_type = 'fixed'

p1 <- expand_parameters(trait_matrix(0.0825, "lma"), p0, FF16_hyperpar,FALSE)
p1$birth_rate <- 20

out <- run_scm_collect(p1)

# select a time slice
i = 120
x <- scm_state(i, out)

test_that("Set patch state", {
  # editable patch object
  types <- extract_RcppR6_template_types(p1, "Parameters")
  patch <- do.call('Patch', types)(p1)

  expect_equal(patch$size, 1) # species
  expect_equal(patch$species[[1]]$size, 0) # cohorts

  patch$introduce_new_cohort(species_index = 1)
  expect_equal(patch$species[[1]]$size, 1) # cohorts

  patch$set_state(x$time, x$species[[1]], n = i, x$env)
  expect_equal(patch$species[[1]]$size, 120)
})

test_that("Set SCM state", {
  types <- extract_RcppR6_template_types(p1, "Parameters")
  scm <- do.call('SCM', types)(p1)
  
  # update introduction schedule
  times <- scm$cohort_schedule$all_times[1]
  expect_equal(length(times[[1]]), 141)
  
  times[[1]] <- c(rep(0, i), times[[1]][-1])
  
  scm$set_cohort_schedule_times(times)
  expect_equal(length(scm$cohort_schedule$all_times[1][[1]]), i + 141 - 1)
  
  # update patch state
  expect_equal(scm$patch$species[[1]]$size, 0)
  
  scm$set_state(0, x$species[[1]], n = i, x$env)
  expect_equal(scm$patch$species[[1]]$size, i)
})


test_that("Set SCM state and collect", {
  y <- run_scm_collect(p1, x)
  
  # check fitness
  expect_equal(y$net_reproduction_ratios, 143.7165, tolerance = 0.0001)

  # check # states, time steps, cohorts
  expect_equal(dim(y$species[[1]]), c(7, 141, 260))
})