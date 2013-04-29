source("helper-tree.R")

context("Misc")

pars <- c(a=1, b=0, c=-5)
x.min <- 0
x.max <- 5
xx <- c(x.min, 2, x.max)
f <- with(as.list(pars), function(x) a*x*x + b*x + c)

## Functor evaluates correctly?
expect_that(tree_module$test_functor(xx, pars),
            equals(f(xx)))

## Find root correctly?
## uniroot(f, c(x.min, x.max))$root

## Analytically
tmp <- with(as.list(pars), (-b + c(-1,1) * sqrt(b^2 - 4*a*c)/(2*a)))
f.root <- tmp[tmp > x.min & tmp < x.max]

## quick sanity check.
expect_that(length(f.root), equals(1))

expect_that(tree_module$test_find_root(pars, x.min, x.max),
            equals(f.root, tolerance=1e-5))

## Root finding to particular value
value <- 10.0
f.value <- uniroot(function(x) f(x) - value, c(x.min, x.max))$root
expect_that(tree_module$test_find_value(pars, value, x.min, x.max),
            equals(f.value, tolerance=1e-5))

## What happens if we're out of range?
expect_that(tree_module$test_find_value(pars, value, x.min, x.min+1),
            throws_error())

## Integration works?
## integrate(f, 0, 5)$value

## Analytically
f.int <- diff(with(as.list(pars), function(x)
                   a/3*x^3 + b/2*x^2 + c*x)(c(x.min, x.max)))

expect_that(tree_module$test_integrator(pars, 0, 5),
            equals(f.int))

set.seed(1)
a <- runif(10)
b <- runif(length(a))
expect_that(tree_module$test_sum_double(a, b), is_identical_to(a + b))

a <- as.integer(a * 10)
b <- as.integer(b * 10)
expect_that(tree_module$test_sum_int(a, b), is_identical_to(a + b))

expect_that(tree_module$test_to_rcpp_matrix(list(a, b)),
            is_identical_to(cbind(a, b, deparse.level=0)))
expect_that(tree_module$test_from_rcpp_matrix(cbind(a, b)),
            is_identical_to(list(a, b)))
