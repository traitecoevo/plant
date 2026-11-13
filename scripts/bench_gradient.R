## Time the reverse sweep of a stand, to compare two builds of plant.
##
## Usage, against an installed plant (one library per build):
##   R_LIBS=<lib> Rscript scripts/bench_gradient.R <model> [reps] [lifetime]
## <model> is FF16 or TF24. Prints one line per repeat: model, repeat, seconds
## for the forward run that records the states, and seconds for the sweep.
##
## Compare builds as scripts/bench_run.R says: alternate invocations, take each
## build's median, and read it beside the instruction count.
##
## The stand is plant's testing parameter set with one species, on the
## birth-date coordinate the sweep requires, differentiated for every metric
## and every column the model declares, through no trade-off function.
library(plant)

args <- commandArgs(trailingOnly = TRUE)
model <- if (length(args) >= 1) args[[1]] else "FF16"
reps <- if (length(args) >= 2) as.integer(args[[2]]) else 5L
lifetime <- if (length(args) >= 3) as.numeric(args[[3]]) else
  c(FF16 = 8, TF24 = 2)[[model]]

p <- default_parameters(model, quiet = TRUE)
p$max_patch_lifetime <- lifetime
p <- add_strategies(p, trait_matrix(0.0825, "lma"), birth_rate = 1.1)
ctrl <- Control(node_density_coordinate = "birth_date")

sweep <- function() {
  t_run <- system.time(
    scm <- run(p, NULL, ctrl, record = "states")$solver)[["elapsed"]]
  t_sweep <- system.time(stand_gradient(scm, hyperpar = NULL))[["elapsed"]]
  c(t_run, t_sweep)
}

invisible(sweep())
for (r in seq_len(reps)) {
  t <- sweep()
  cat(sprintf("%s\t%d\t%.4f\t%.4f\n", model, r, t[[1]], t[[2]]))
}
