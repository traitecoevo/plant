test_that("plant_log_console is deprecated and does nothing", {
  expect_warning(out <- plant_log_console(), "no effect")
  expect_null(out)
})
