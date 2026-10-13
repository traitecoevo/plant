## Capture the outputs of runs that use none of the opt-in features, for
## tests/testthat/test-default-baseline.R. Those runs must not move when an
## opt-in feature is added, so the baseline is captured once, from a build that
## predates the features it guards, and is re-blessed only for a deliberate
## change to the default path.
##
## Usage, against an installed plant:
##   Rscript scripts/capture-default-baseline.R [out.rds]
## The cases are built by default_baseline_cases(), defined in the test's
## helper, so the test and the capture run the same code. It uses only API
## that predates #661, so it runs against the commit before it (9f04309b),
## which is where the recorded baseline comes from.

args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args) > 0) args[[1]] else
  "tests/testthat/test_data/default_baseline.rds"
library(plant)
message("plant from ", find.package("plant"))
source("tests/testthat/helper-default-baseline.R")
cases <- default_baseline_cases()
attr(cases, "plant") <- as.character(packageVersion("plant"))
attr(cases, "platform") <- paste(R.version$platform, R.version$version.string)
saveRDS(cases, out)
message("wrote ", length(cases), " cases to ", out)
