## Time one stand run, to compare two builds of plant.
##
## Usage, against an installed plant (one library per build):
##   R_LIBS=<lib> Rscript scripts/bench_run.R <model> [reps] [lifetime]
## <model> is FF16 or TF24. Prints one line per repeat: model, repeat, seconds.
## A repeat is ten runs for FF16, whose single run is too short to time, and
## one for TF24.
##
## Two builds cannot load into one R session, so compare them by alternating
## invocations, base then head, several rounds, and take each build's median.
## Wall-clock on a laptop swings by more than most changes are worth, so read
## the medians beside the instruction count `/usr/bin/time -l` reports for the
## same invocation (see .claude/measurements.md).
##
## The stand is plant's testing parameter set with one species, run on the
## default schedule with no refinement, so the work is the same in both builds.
library(plant)

args <- commandArgs(trailingOnly = TRUE)
model <- if (length(args) >= 1) args[[1]] else "FF16"
reps <- if (length(args) >= 2) as.integer(args[[2]]) else 5L
lifetime <- if (length(args) >= 3) as.numeric(args[[3]]) else
  c(FF16 = 105, TF24 = 10)[[model]]

p <- default_parameters(model, quiet = TRUE)
p$max_patch_lifetime <- lifetime
p <- add_strategies(p, trait_matrix(0.0825, "lma"), birth_rate = 20)
env <- Environment(model)
ctrl <- Control()

runs <- c(FF16 = 10L, TF24 = 1L)[[model]]
invisible(run(p, env, ctrl))
for (r in seq_len(reps)) {
  t <- system.time(for (i in seq_len(runs)) run(p, env, ctrl)$solver)[["elapsed"]]
  cat(sprintf("%s\t%d\t%.4f\n", model, r, t))
}
