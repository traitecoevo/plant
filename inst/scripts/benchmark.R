#!/usr/local/bin/Rscript

# Run benchmarks using package bench
library(bench)

# load plant package
devtools::load_all()

bench_scm <- function(x) {

	p0 <- scm_base_parameters(x)
	p <- expand_parameters(trait_matrix(0.0825, "lma"), p0, mutant=FALSE)
	if(grepl("K93", x))
      p$k_I <- 1e-3
	res <- run_scm(p)
}

results <- 
	bench::press(
  x = c("FF16", "FF16r", "K93"),
  {
    bench::mark(
      min_iterations = 2,
      bench_scm(x)
    )
  }
)

results
