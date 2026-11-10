# What every model's individual must provide, checked once over all of them.

models <- c("FF16", "K93", "TF24", "TF24f", "TF24t", "TF24w")

test_that("the first three states are height, mortality and fecundity", {
  # The core reads these three by fixed index (HEIGHT_INDEX and the rest), so
  # every model must declare them first and in this order.
  for (m in models) {
    ind <- get(paste0(m, "_Individual"))()
    expect_identical(ind$ode_names[1:3], c("height", "mortality", "fecundity"),
                     info = m)
  }
})

test_that("an individual holds its strategy, and reads height by name and by slot", {
  for (m in models) {
    s <- get(paste0(m, "_Strategy"))()
    ind <- get(paste0(m, "_Individual"))(s)
    expect_identical(ind$strategy, s, info = m)
    ind$set_state("height", 10)
    expect_identical(ind$state("height"), 10, info = m)
    expect_identical(ind$state("height"),
                     ind$internals$states[which(ind$ode_names == "height")],
                     info = m)
  }
})
