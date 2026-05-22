
context("Extrinsic drivers")

test_that("ExtrinsicDrivers can be instantiated", {
  
  expect_silent(drv <- ExtrinsicDrivers())

  expect_contains(names(drv), c('clone', 'evaluate', 'evaluate_range', 'get_names', 'initialize', 'set_constant', 'set_extrapolate', 'set_variable'))
  
  expect_silent(drv$set_constant("test", 1000))
 
  # extract
  expect_equal(drv$evaluate("test", 1), 1000)
  expect_equal(drv$evaluate("test", 10), 1000)

  # set new value
  expect_silent(drv$set_constant("test", 1.5))
  expect_equal(drv$evaluate("test", 1), 1.5)
  expect_equal(drv$evaluate("test", 10), 1.5)

  # set a bunch of values
  for(i in 1:10)
    expect_silent(drv$set_constant(letters[i], i))
  for (i in 1:10)
    expect_equal(drv$evaluate(letters[i], 1), i)

  # set new value
  for(i in 1:10)
    expect_silent(drv$set_constant(letters[i], i-1))
  for (i in 1:10)
    expect_equal(drv$evaluate(letters[i], 1), i-1)
 
  # expected erors
  expect_error(drv$evaluate("wrong", 1))
  
  # a function
  x <- 1:10
  y <- sin(x)
  expect_silent(drv$set_variable("sin", x,y))
  for (i in 1:10) 
    expect_equal(drv$evaluate("sin", i), y[i])
})
