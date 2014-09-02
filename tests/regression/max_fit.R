library(tree)
library(testthat)

p <- ebt_base_parameters()
r <- tree:::max_growth_rate("lma", 0.2, p)

# Compare against default value, claulcated with version db85b05c9236836a640c99a6820a22a5da8603c1
expect_that(r,equals(11.16109, tolerance=1e-4))
