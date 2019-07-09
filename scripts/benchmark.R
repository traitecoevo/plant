
# Run benchmarks using package microbenchmark
library(microbenchmark)

# load plant package
pkgload::load_all()

# set parameters
p0 <- scm_base_parameters("FF16")
p <- expand_parameters(trait_matrix(0.0825, "lma"), p0, FALSE)
p$seed_rain <- 20 # close to equilibrium

# benchmark using microbenchmark
res <- microbenchmark(run_scm(p), times=5L, unit="s")

# Print results
print(res)

#  ggplot2::autoplot(res)

