## Writes inst/scenarios/scenario_rainfall_drivers.rds: the seasonal rainfall
## series the TF24 scenario gateway runs on, one table per (rainfall_mean,
## rainfall_amp_frac) pair the scenario mapping uses.
##
## This is the sinusoid build_scenario() used to generate inline, kept as a
## fixture so the gateway baseline (tests/testthat/test_data/scenario_baseline.rds)
## stays reproducible after the generator moved to fieldtrial. The arithmetic
## below is the former code verbatim, so the knots are bit-identical to what
## the baseline was recorded with: mean * (1 + amp * sin(2 pi t)), trough
## clamped at zero, 48 knots per year (at 6 per year a cubic through the sine
## was 0.18% off at the peak; at 48 the error is ~1e-6).
##
## Run from the package root: Rscript scripts/make-scenario-drivers.R
## Re-run only when the mapping gains a new rainfall level or amplitude, and
## say so in the PR: a changed fixture is a changed gateway.

max_patch_lifetime <- 100
points_per_year <- 48L
mapping <- utils::read.csv("inst/scenarios/scenario_mapping.csv",
                           stringsAsFactors = FALSE, check.names = FALSE)
means <- as.numeric(mapping$value[mapping$param == "rainfall_mean"])
amps <- as.numeric(mapping$value[mapping$param == "rainfall_amp_frac"])
amps <- amps[amps != 0]

t <- seq(0, max_patch_lifetime,
         length.out = max(200L, as.integer(ceiling(max_patch_lifetime *
                                                     points_per_year))))
out <- list()
for (mean_rain in means) {
  for (amp_frac in amps) {
    y <- mean_rain * (1 + amp_frac * sin(2 * pi * t))
    y <- pmax(y, 0)
    out[[sprintf("mean=%s;amp=%s", format(mean_rain), format(amp_frac))]] <-
      data.frame(time = t, rainfall = y)
  }
}
attr(out, "max_patch_lifetime") <- max_patch_lifetime
attr(out, "points_per_year") <- points_per_year
saveRDS(out, "inst/scenarios/scenario_rainfall_drivers.rds", version = 3)
cat("wrote", length(out), "series:", paste(names(out), collapse = ", "), "\n")
