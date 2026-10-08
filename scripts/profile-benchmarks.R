#!/usr/bin/env Rscript

# Profile the same SCM scenarios measured by run_plant_benchmarks().
#
# Usage:
#   Rscript scripts/profile-benchmarks.R
#   Rscript scripts/profile-benchmarks.R FF16 K93
#
# Run `make compile` first; this script loads the already-built package DLL
# rather than recompiling through load_all().
#
# Outputs are written below tmp/profile-benchmarks/<timestamp>/ so they do not
# clutter version control.

args <- commandArgs(trailingOnly = TRUE)
strategies <- if (length(args)) args else "FF16"
iterations <- as.integer(Sys.getenv("PLANT_BENCHMARK_ITERATIONS", "1"))
rprof_interval <- as.numeric(Sys.getenv("PLANT_RPROF_INTERVAL", "0.005"))
profile_repeats <- as.integer(Sys.getenv("PLANT_PROFILE_REPEATS", "1"))
sample_seconds <- as.integer(Sys.getenv("PLANT_SAMPLE_SECONDS", "0"))

timestamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
out_dir <- file.path("tmp", "profile-benchmarks", timestamp)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Loading plant from existing compiled DLL")
pkgload::load_all(compile = FALSE, debug = FALSE, quiet = TRUE)

strategy_constructors <- list(
  FF16 = FF16_Strategy,
  TF24 = TF24_Strategy,
  K93 = K93_Strategy
)

unknown <- setdiff(strategies, names(strategy_constructors))
if (length(unknown)) {
  stop("Unknown strategy type(s): ", paste(unknown, collapse = ", "), call. = FALSE)
}

profile_cases <- list(
  scm = function(strategy) {
    p0 <- scm_base_parameters(strategy)
    p <- add_strategies(p0, trait_matrix(0.0825, "lma"))
    run_scm(p)
    invisible(NULL)
  },
  build_schedule = function(strategy) {
    p <- scm_base_parameters(strategy)
    p$strategies <- list(strategy_constructors[[strategy]]())
    p$birth_rate <- 0.1
    run_scm(p, refine_schedule = TRUE)
    invisible(NULL)
  }
)

write_benchmark_summary <- function(x, path) {
  summary <- data.frame(
    expression = as.character(x$expression),
    strategy = x$strategy,
    min = as.character(x$min),
    median = as.character(x$median),
    itr_sec = x$`itr/sec`,
    mem_alloc = as.character(x$mem_alloc),
    n_itr = x$n_itr,
    n_gc = x$n_gc,
    total_time = as.character(x$total_time)
  )
  utils::write.csv(summary, path, row.names = FALSE)
}

write_rprof_summary <- function(profile_file, prefix) {
  summary <- summaryRprof(profile_file)

  by_self <- as.data.frame(summary$by.self)
  by_total <- as.data.frame(summary$by.total)

  if (nrow(by_self)) {
    by_self$function_name <- rownames(by_self)
    by_self <- by_self[, c("function_name", setdiff(names(by_self), "function_name"))]
  }
  if (nrow(by_total)) {
    by_total$function_name <- rownames(by_total)
    by_total <- by_total[, c("function_name", setdiff(names(by_total), "function_name"))]
  }

  utils::write.csv(by_self, paste0(prefix, "-by-self.csv"), row.names = FALSE)
  utils::write.csv(by_total, paste0(prefix, "-by-total.csv"), row.names = FALSE)
}

run_with_optional_sample <- function(expr, sample_file) {
  if (sample_seconds <= 0) {
    force(expr)
    return(invisible(NULL))
  }

  if (!file.exists("/usr/bin/sample")) {
    warning("PLANT_SAMPLE_SECONDS was set, but /usr/bin/sample was not found")
    force(expr)
    return(invisible(NULL))
  }

  parent_pid <- Sys.getpid()
  sample_log <- sub("[.]sample[.]txt$", ".sample.log", sample_file)
  sampler <- parallel::mcparallel({
    output <- system2(
      "/usr/bin/sample",
      c(as.character(parent_pid), as.character(sample_seconds),
        "-file", sample_file),
      stdout = TRUE,
      stderr = TRUE
    )
    status <- attr(output, "status")
    if (is.null(status)) {
      status <- 0L
    }
    writeLines(output, con = sample_log, useBytes = TRUE)
    if (!identical(status, 0L)) {
      cat("sample exited with status", status, "\n",
          file = sample_log, append = TRUE)
    }
    status
  })

  force(expr)
  parallel::mccollect(sampler)
  invisible(NULL)
}

run_profile_case <- function(case, strategy) {
  for (i in seq_len(profile_repeats)) {
    case(strategy)
  }
  invisible(NULL)
}

message("Writing outputs to ", out_dir)

benchmarks <- run_plant_benchmarks(
  strategy_types = strategy_constructors[strategies],
  iterations = iterations
)
write_benchmark_summary(benchmarks, file.path(out_dir, "benchmark-summary.csv"))
saveRDS(benchmarks, file.path(out_dir, "benchmark-summary.rds"))

for (strategy in strategies) {
  for (case_name in names(profile_cases)) {
    prefix <- file.path(out_dir, paste(strategy, case_name, sep = "-"))
    profile_file <- paste0(prefix, ".Rprof")

    message("Profiling ", strategy, " / ", case_name)
    gc()
    sample_file <- paste0(prefix, ".sample.txt")
    Rprof(profile_file, interval = rprof_interval, memory.profiling = TRUE)
    run_with_optional_sample(
      run_profile_case(profile_cases[[case_name]], strategy),
      sample_file
    )
    Rprof(NULL)

    write_rprof_summary(profile_file, prefix)
  }
}

message("Done")
