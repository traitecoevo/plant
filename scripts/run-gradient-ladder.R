# Run the checks on the reverse-mode gradient, and print where they stand.
#
#   Rscript scripts/run-gradient-ladder.R                 every file
#   Rscript scripts/run-gradient-ladder.R one-cohort       one file
#
# Run from the package root, and build optimised first: an unoptimised build
# makes the trajectory fixtures several times slower.
#
#   make
#
# ---------------------------------------------------------------------------
#
# WHAT THESE CHECKS ARE, AND WHY THERE ARE SO MANY OF THEM
#
# The claim under test is that the reverse sweep is the exact transpose of the
# forward run. That claim cannot be checked by a finite difference, because
# differencing is the thing the sweep replaces -- and near a coincidence a
# difference does not merely lose accuracy, it fails to converge. So the claim is
# checked against references that share no code with the sweep.
#
# ⚠️ THE REFERENCES ARE ENUMERATED IN `docs/design/reverse-mode.md` AND NOT HERE.
# R1 a forward tangent, R2 the whole Jacobian at one cohort, R3 a difference of
# whole runs, R4 the model rebuilt from its parameters -- each with what it is
# blind to, which is the half that decides why four are needed. This file maps
# files to those numbers and restates none of them: the two lists were written
# separately once and had drifted into two different sets of four.
#
# The fixtures form a progression -- the "ladder" the filenames name -- from the
# smallest stand whose Jacobian can be formed whole, through accumulation across
# cohorts and species, to a stand whose state vector changes width. Each step adds
# one mechanism and uses the strongest reference its size still permits:
#
#   no reference needed   model-invariants, identity
#   R2                    one-cohort, switches
#   R1                    two-species, columns, introductions, first-range, recruit
#   R4                    factorisation
#   R3                    whole-run-difference, declared-zero
#   whether they bite     sweep, injection
#
# ⚠️ A MARGIN IS NOT A PASS. Every check reports how much of its tolerance budget
# it used, because a check passing at three-quarters of budget is a check about to
# stop working. And a suite of margins says nothing about whether the checks would
# notice a defect at all, which is what the fault injections establish: two of the
# first three injections tried here failed to fail.
#
# ⚠️ A SKIP IS NOT A PASS EITHER, and this suite has none. Every gate that could
# produce one has been removed, so a skip appearing here is a check that has
# stopped asking its question rather than a condition it declined to meet. Read
# the skip count beside the failures every time.

suppressMessages({
  library(odelia)
  pkgload::load_all(".", quiet = TRUE)
  library(testthat)
})

# The order a reader should meet them: smallest fixture whose Jacobian can be
# formed whole, then accumulation, then a state vector that changes width. ONLY
# the order lives here. What each file claims is the first line of that file, read
# below -- restating it here made two sources that had to be kept equal by hand,
# and they had drifted into two different sentences.
ladder <- c("model-invariants", "identity",
            "one-cohort", "switches",
            "two-species", "columns", "introductions", "first-range", "recruit",
            "factorisation",
            "whole-run-difference", "declared-zero",
            "sweep", "injection")

# ⚠️ CHECKED AGAINST THE DIRECTORY, because a list read off a directory loses the
# order and a list written by hand goes silently incomplete when the directory
# grows. A file missing from `ladder` would never run under this script and would
# look like a suite that passes.
claim_of <- function(name, dir) {
  first <- readLines(file.path(dir, paste0("test-gradient-ladder-", name, ".R")),
                     n = 1L, warn = FALSE)
  sub("\\.$", "", sub("^# *", "", first))
}

dir <- file.path("tests", "testthat")

# ⚠️ MAINTAINED, AND CHECKED, because a list read off a directory goes silently
# incomplete when the directory grows. A file with no entry here would otherwise
# never run under this script and would look like a suite that passes.
on_disk <- sub("^test-gradient-ladder-", "",
               sub("\\.R$", "",
                   basename(Sys.glob(file.path(dir, "test-gradient-ladder-*.R")))))
undescribed <- setdiff(on_disk, ladder)
if (length(undescribed) > 0L) {
  stop("These files are not in `ladder`, so this script would not run them: ",
       paste(undescribed, collapse = ", "), ". Add each where a reader should ",
       "meet it.")
}
missing <- setdiff(ladder, on_disk)
if (length(missing) > 0L) {
  stop("These entries in `ladder` name no file: ",
       paste(missing, collapse = ", "))
}

selected <- commandArgs(trailingOnly = TRUE)
if (length(selected) == 0L) selected <- ladder
unknown <- setdiff(selected, ladder)
if (length(unknown) > 0L) {
  stop("No such file: ", paste(unknown, collapse = ", "),
       ". Available: ", paste(ladder, collapse = ", "))
}

results <- list()
for (name in selected) {
  path <- file.path(dir, paste0("test-gradient-ladder-", name, ".R"))
  cat(sprintf("\n===== %s -- %s\n", name, claim_of(name, dir)))
  out <- testthat::test_file(path, package = "plant", reporter = "summary")
  frame <- as.data.frame(out)
  results[[name]] <- c(
    pass = sum(frame$passed),
    fail = sum(frame$failed),
    skip = sum(frame$skipped),
    error = sum(frame$error))
}

cat("\n\n===== where the checks stand =====\n")
table <- do.call(rbind, results)
print(table)

blocked <- rownames(table)[table[, "skip"] > 0]
if (length(blocked) > 0) {
  cat("\nSkips name something these checks cannot yet ask, and they are not\n",
      "passes. Skipped: ", paste(blocked, collapse = ", "), "\n", sep = "")
}
