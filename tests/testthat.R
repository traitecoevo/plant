library(testthat)
library(plant)

# PLANT_TEST_EXCLUDE, when set, is a regular expression matched against each
# test file's name without "test-" and ".R"; matching files are not run. CI
# sets it on the Windows leg of a pull request to leave out the slowest files.
exclude <- Sys.getenv("PLANT_TEST_EXCLUDE")
if (nzchar(exclude)) {
  message("Leaving out test files matching ", exclude)
  test_check("plant", filter = exclude, invert = TRUE)
} else {
  test_check("plant")
}
