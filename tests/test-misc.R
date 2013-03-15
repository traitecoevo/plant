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

## Integration works?
## integrate(f, 0, 5)$value

## Analytically
f.int <- diff(with(as.list(pars), function(x)
                   a/3*x^3 + b/2*x^2 + c*x)(c(x.min, x.max)))

expect_that(tree_module$test_integrator(pars, 0, 5),
            equals(f.int))
