# The differenced reference the analytic leaf rows are refereed against.
#
# WHY THIS EXISTS AND WHY IT COMES FIRST. Report 10 replaces differenced leaf rows
# with analytic ones. A change from differenced to analytic returns a finite,
# plausible, wrong gradient when it is wrong, and the differenced implementation is
# the only instrument that shares no assumption with the sweep -- so it is the
# reference, and it stops existing the moment it is deleted. This captures it while
# it is still there.
#
# WHAT IT IS. `ladder_run_difference_pair` rebuilds the strategy and re-runs a
# stand where two species COMPETE, so the row it reads carries the belowground
# coupling a single-species reference cannot see: rows six per cent out have agreed
# with a one-species reference to 1e-06.
#
# ⚠️ THE STEP IS NOT A DETAIL. A column whose sensitivity is an order below its
# neighbour's sits close to this difference's own floor, and its reading then moves
# with the step instead of holding. Every step's reading is written out, not just
# the chosen one, so a later reader can see whether a column had converged rather
# than trusting that it had.
#
# ⚠️ AND THE REGIME DECIDES WHICH KIND OF OPERATING POINT IS EXERCISED. A reference
# taken on one regime says nothing about the pinned or shut branches. Which kinds
# each regime actually reached is recorded beside the numbers, by the model's own
# counter, rather than assumed from the regime's name.
#
#   [REFERENCE_WORKERS=n] Rscript scripts/capture-reference-gradient.R [outdir]
#
# Run from the package root. Writes <outdir>/reference-gradient.tsv and
# <outdir>/reference-kinds.tsv, both %.17g so a double round-trips; the default
# outdir is where the whole-run rung reads them from.

root <- normalizePath(".")
outdir <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[[1]] else
  "tests/testthat/reference"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
outdir <- normalizePath(outdir)

library(odelia)
pkgload::load_all(root, quiet = TRUE)
# Sourced from inside the test directory, because one helper resolves the
# testthat root when it is read and finds nothing from anywhere else. Paths above
# are already absolute for the same reason.
setwd(file.path(root, "tests", "testthat"))
for (f in list.files(".", "^helper.*\\.R$")) {
  source(f)
}
steps <- c(1e-6, 1e-5, 1e-4, 1e-3)
workers <- as.integer(Sys.getenv("REFERENCE_WORKERS",
                                 unset = max(1L, parallel::detectCores() - 1L)))
times <- list(c(0, 0.63), c(0, 0.41))
regimes <- ladder_reference_regimes()

g17 <- function(x) formatC(x, format = "g", digits = 17)

# ⚠️ THE SAME EXPRESSION THE WHOLE-RUN RUNG ASSERTS THIS FILE COVERS, so the
# capture and the check cannot disagree about what a column is. Written two ways
# they did: the rung read the strategy's registered parameters and this read the
# names a patch carried, and seven columns sat uncaptured and unrefereed.
columns <- setdiff(names(TF24_Strategy()$pars),
                   names(census_undifferentiable_tf24()))
metrics <- census_metric_names_tf24()
cat(sprintf("columns %d   metrics %d   regimes %d   species 2\n",
            length(columns), length(metrics), length(regimes)))

kinds <- list()
rows <- list()

for (rg in regimes) {
  # One run per regime, before any differencing, to record what it reached. The
  # counter is cleared first: it accumulates over every solve the process has done.
  p <- ladder_parameters(c("fast", "slow"), k_I = rg$k_I)
  p$node_schedule_times <- times
  scm <- ladder_run(p, env = ladder_environment(rg$rain,
                                                if (is.null(rg$amplitude)) 0
                                                else rg$amplitude))
  counts <- census_operating_point_counts_tf24(scm)
  for (s in seq_along(counts)) {
    kinds[[length(kinds) + 1L]] <- data.frame(
      regime = rg$name, species = s,
      kind = census_operating_point_names_tf24(), count = counts[[s]],
      stringsAsFactors = FALSE)
  }
  reached <- census_operating_point_names_tf24()[
    Reduce(`+`, counts) > 0]
  cat(sprintf("%-9s reached: %s\n", rg$name, paste(reached, collapse = " ")))
  flush(stdout())

  for (sp in 1:2) {
    # A column the strategy does not carry has nothing to difference, and a
    # column whose difference cannot be taken is recorded as such rather than
    # dropped -- an absent row and a refused one are not the same statement.
    one <- function(nm) {
      got <- tryCatch(
        ladder_run_difference_pair(nm, species = sp, steps = steps,
                                   times = times, regime = rg),
        error = function(e) conditionMessage(e))
      head <- data.frame(regime = rg$name, species = sp, parameter = nm,
                         stringsAsFactors = FALSE)
      if (is.character(got)) {
        return(data.frame(head, metric = NA_character_, s1 = NA, s2 = NA,
                          s3 = NA, s4 = NA, converged = NA, step = NA,
                          spread = NA, note = got, stringsAsFactors = FALSE))
      }
      data.frame(head, metric = rownames(got$values), s1 = got$values[, 1],
                 s2 = got$values[, 2], s3 = got$values[, 3],
                 s4 = got$values[, 4], converged = unname(got$gradient),
                 step = got$step, spread = got$spread, note = "",
                 stringsAsFactors = FALSE)
    }
    # The columns are independent and there are 480 of them at nine runs each, so
    # they go out to separate processes -- but a forked worker that dies reports
    # only that all cores errored, so each carries its own condition back.
    got <- parallel::mclapply(columns,
                              function(nm) tryCatch(one(nm), condition = identity),
                              mc.cores = workers)
    for (i in seq_along(got)) if (inherits(got[[i]], "condition")) {
      stop(sprintf("%s / %d.%s: %s", rg$name, sp, columns[[i]],
                   conditionMessage(got[[i]])), call. = FALSE)
    }
    rows <- c(rows, got)
    cat(sprintf("%-9s species %d done\n", rg$name, sp)); flush(stdout())
  }
}

out <- do.call(rbind, rows)
num <- c("s1", "s2", "s3", "s4", "converged", "step", "spread")
for (v in num) out[[v]] <- g17(out[[v]])
write.table(out, file.path(outdir, "reference-gradient.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
write.table(do.call(rbind, kinds), file.path(outdir, "reference-kinds.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

refused <- sum(nzchar(out$note))
cat(sprintf("\nwrote %d rows (%d could not be differenced) to %s\n",
            nrow(out), refused, outdir))
# The coverage statement, by name. A kind no regime reached is what a later
# reader needs to know before trusting the reference on that branch.
all_kinds <- census_operating_point_names_tf24()
seen <- unique(do.call(rbind, kinds)[do.call(rbind, kinds)$count > 0, "kind"])
cat("kinds reached: ", paste(seen, collapse = " "), "\n")
cat("kinds NOT reached: ",
    paste(setdiff(all_kinds, seen), collapse = " "), "\n")
