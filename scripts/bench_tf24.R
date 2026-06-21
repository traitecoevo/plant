# Reproducible cumulative benchmark for TF24 root-water-uptake optimisations.
# Loads the already-make-built .so WITHOUT recompiling (so we measure the
# compiled binary, not a load_all() rebuild — see the build note in
# test-root.qmd), times a single-species MPL=10 SCM run, and reports a numeric
# trajectory metric. Build both sides clean with `make clean && make compile`
# first; set BENCH_REPS to change the number of timed reps (default 5).

suppressMessages(pkgload::load_all(".", recompile = FALSE, quiet = TRUE))

p0 <- scm_base_parameters("TF24")
p0$max_patch_lifetime <- 10
env  <- Environment("TF24")
ctrl <- control()
p1 <- expand_parameters(trait_matrix(0.0825, "lma"), p0, TF24_hyperpar,
                        birth_rate_list = list(20))

# warm-up (JIT/allocator/caches)
invisible(run_scm(p1, env, ctrl))

reps <- as.integer(Sys.getenv("BENCH_REPS", "5"))
t <- numeric(reps)
op <- NA_real_
for (i in seq_len(reps)) {
  tt <- system.time({ out <- run_scm(p1, env, ctrl) })
  t[i] <- tt[["elapsed"]]
  op <- out$offspring_production
}

cat(sprintf("REPS=%d\n", reps))
cat(sprintf("elapsed_s: min=%.3f median=%.3f mean=%.3f max=%.3f\n",
            min(t), median(t), mean(t), max(t)))
cat(sprintf("offspring_production: %.10g\n", op))
cat("raw:", paste(sprintf("%.3f", t), collapse=" "), "\n")
