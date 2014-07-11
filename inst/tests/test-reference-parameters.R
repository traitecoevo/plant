source("helper-tree.R")

skip <- file.exists(".SKIP_REFERENCE_TESTS")

## Disabling these tests for now: we've diverged in physiological
## model and in the number of parameters.  I will support putting
## these back together at some point, but for now leaving things be.
skip <- TRUE

## This is basically a test of if we're running in Rich's computer.
if (!skip) {

context("Reference Model Parameters")

## Path with some example parameter files:
path.src <- system.file("reference/falster-traitdiversity",
                        package="tree")

## Path to store output:
path.cmp <- tempfile("cmp")

## 1: Try to read in the parameters:
obj <- tree:::read.reference.parameters(path.src)

## 2: Try to write the parameters back out
tree:::write.reference.parameters(obj, path.cmp)

## Compare these:
obj.cmp <- tree:::read.reference.parameters(path.cmp)

test_that("Read/write of reference objects doesn't alter them", {
  expect_that(obj.cmp, equals(obj))
})

## 3: Try pulling into a tree object:
obj.tree <- tree:::parameters.from.reference(obj)
test_that("Conversion to tree object doesn't alter parameters", {
  expect_that(tree:::reference.from.parameters(obj.tree),
              equals(obj))
})

## Public interface:
obj <- load.reference.parameters(path.src)
save.reference.parameters(obj, path.cmp)
obj.cmp <- load.reference.parameters(path.cmp)

test_that("Complete save/load doesn't alter parameters", {
  expect_that(tree:::reference.compare.parameters(obj.cmp, obj),
              is_true())
})

## Now, try building some objects and saving them.  First, let's build
## an object that cannot be losslessly converted to run with the
## reference model:
p <- new(Parameters)
p$add_strategy(new(Strategy, list(rho=0.1)))
p$add_strategy(new(Strategy, list(rho=0.3)))
p$seed_rain <- c(.1, .2)
test_that("Can't create object with variable rho", {
  expect_that(save.reference.parameters(p, path.cmp),
              throws_error())
})

## With a normal parameter:
p <- new(Parameters)
p$add_strategy(new(Strategy, list(a1=0.1)))
p$add_strategy(new(Strategy, list(a1=0.3)))
p$seed_rain <- c(.1, .2)
test_that("Can't create object with variable parameter (a1)", {
  expect_that(save.reference.parameters(p, path.cmp),
              throws_error())
})

## And then build one that we *can* tweak:
p <- new(Parameters)
p$add_strategy(new(Strategy, list(lma=0.1)))
p$add_strategy(new(Strategy, list(lma=0.3)))
p$seed_rain <- c(.1, .2)
save.reference.parameters(p, path.cmp)

p.cmp <- load.reference.parameters(path.cmp)
test_that("Complete save/load doesn't alter parameters", {
  expect_that(tree:::reference.compare.parameters(p.cmp, p),
              is_true())
})

## This is not run by default because it's set up for my machine only.
## Still thinking about portability and long term comparison.
if (FALSE) {
  path.evolve <- "../../../falster-traitdiversity/src"
  run.reference(path.cmp, p, path.evolve=path.evolve, verbose=FALSE)

  output <- load.reference.output(path.cmp)
  test_that("Output contains correct parameters", {
    expect_that(tree:::reference.compare.parameters(output$parameters, p),
                is_true())
  })

  ## This is a bit mean, but checks that the output is what we expect.
  library(digest)
  test_that("Output is as expected (harsh test)", {
    expect_that(digest(output[names(output) != "parameters"], "sha1"),
                is_identical_to("725c32f87efa820f031e7c56f9564f024bc1d876"))
  })
}
}
