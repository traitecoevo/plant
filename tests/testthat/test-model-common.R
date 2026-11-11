# What every model's individual must provide, checked once over all of them.

test_that("the first three states are height, mortality and fecundity", {
  # The core reads these three by fixed index (HEIGHT_INDEX and the rest), so
  # every model must declare them first and in this order.
  for (m in all_models) {
    ind <- get(paste0(m, "_Individual"))()
    expect_identical(ind$ode_names[1:3], c("height", "mortality", "fecundity"),
                     info = m)
  }
})

test_that("an individual holds its strategy, and reads height by name and by slot", {
  for (m in all_models) {
    s <- get(paste0(m, "_Model"))()
    ind <- get(paste0(m, "_Individual"))(s)
    expect_identical(ind$strategy, s, info = m)
    ind$set_state("height", 10)
    expect_identical(ind$state("height"), 10, info = m)
    expect_identical(ind$state("height"),
                     ind$internals$states[which(ind$ode_names == "height")],
                     info = m)
  }
})
