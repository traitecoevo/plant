## Capture the outputs of runs that use none of the opt-in features, for
## tests/testthat/test-default-baseline.R. Those runs must not move when an
## opt-in feature is added, so the baseline is captured once, from a build that
## predates the features it guards, and is re-blessed only for a deliberate
## change to the default path.
##
## Usage, against an installed plant:
##   Rscript scripts/capture-default-baseline.R [out.rds]
## The cases are built by default_baseline_cases(), defined in the test's
## helper, so the test and the capture run the same code. The helper runs
## through run(), so the capture needs a build that has it (#428 onwards). The
## recorded baseline comes from the commit before #661 (9f04309b), captured
## when the helper still used the API of that commit; to re-capture from a
## build that old, use the helper as it stood at that commit.

args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args) > 0) args[[1]] else
  "tests/testthat/test_data/default_baseline.rds"
library(plant)
message("plant from ", find.package("plant"))
## helper-plant.R defines testing_parameters() and ff16_parameters(), which
## the cases use.
source("tests/testthat/helper-plant.R")
source("tests/testthat/helper-default-baseline.R")
cases <- default_baseline_cases()
attr(cases, "plant") <- as.character(packageVersion("plant"))
attr(cases, "platform") <- paste(R.version$platform, R.version$version.string)
saveRDS(cases, out)
message("wrote ", length(cases), " cases to ", out)
