source("helper-tree.R")

context("Misc")

pars <- c(1, 0, -5)
expect_that(tree_module$test_functor(c(0, 2, 5), pars),
            equals(c(-5, -1, 20)))

expect_that(tree_module$test_find_root(pars, 0, 5),
            equals(sqrt(5), tolerance=1e-5))
