# Snapshot / check the AD-gradient correctness-regression fixture (#472 scope B,
# refactor+optimize phase). The "AD-vs-AD baseline": pins each validated Jacobian
# to a committed snapshot so the refactor (harvest -> C++, engine unification) is
# provably value-preserving. See notes/ad-refactor-optimize-roadmap.md and
# tests/testthat/helper-gradient-fixture.R (the spec, shared with the testthat test).
#
# Usage (build optimized first: make compile):
#   Rscript --no-init-file scripts/gradient_fixture.R snapshot   # write the .rds (PRE-refactor)
#   Rscript --no-init-file scripts/gradient_fixture.R check      # compare current tree to the .rds
#
# `snapshot` is run ONCE on the pre-refactor HEAD and the .rds committed; `check`
# is the same comparison the testthat test runs, handy for a quick ad-hoc gate.

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[[1]] else "check"
path <- if (length(args) > 1) args[[2]] else "."
suppressMessages(pkgload::load_all(path, compile = FALSE, quiet = TRUE))
source(file.path(path, "tests/testthat/helper-gradient-fixture.R"))

rds <- file.path(path, "tests/testthat/fixtures/gradient-baseline.rds")
specs <- gradient_fixture_specs()

if (identical(mode, "snapshot")) {
  dir.create(dirname(rds), showWarnings = FALSE, recursive = TRUE)
  baseline <- lapply(specs, function(s) gradient_fixture_flatten(s$compute()))
  saveRDS(baseline, rds, version = 2)
  cat(sprintf("snapshot: wrote %d fixtures to %s\n", length(baseline), rds))
  for (nm in names(baseline))
    cat(sprintf("  %-20s n=%d  sum|.|=%.10g\n", nm, length(baseline[[nm]]),
                sum(abs(baseline[[nm]]))))
} else {
  if (!file.exists(rds)) stop("no baseline at ", rds, " -- run snapshot first")
  baseline <- readRDS(rds)
  fail <- 0L
  for (nm in names(specs)) {
    tier <- specs[[nm]]$tier
    tol  <- gradient_fixture_tol(tier)
    cur  <- gradient_fixture_flatten(specs[[nm]]$compute())
    ref  <- baseline[[nm]]
    if (is.null(ref)) { cat(sprintf("%-20s  NEW (no baseline entry)\n", nm)); next }
    rel  <- abs(cur - ref) / pmax(abs(ref), 1e-300)
    rel[!is.finite(rel)] <- ifelse(cur[!is.finite(rel)] == ref[!is.finite(rel)], 0, Inf)
    worst <- max(rel)
    ok <- worst <= tol
    if (!ok) fail <- fail + 1L
    cat(sprintf("%-20s  tier=%-5s  max_rel=%.3e  tol=%.0e  %s\n",
                nm, tier, worst, tol, if (ok) "PASS" else "FAIL"))
  }
  if (fail > 0L) { cat(sprintf("\n%d fixture(s) FAILED\n", fail)); quit(status = 1L) }
  cat("\nall fixtures PASS\n")
}
