# Same-machine A/B timing harness for FF16 SCM, used to compare two builds of
# plant (e.g. this branch vs another branch, or before/after a change).
#
# You cannot load two builds of `plant` in one R session, so build each and run
# this script once per build in a separate Rscript, interleaving the runs to
# control for machine drift. Reports the one-iteration median (ms, 7 reps) and
# the 20-repeat wall-clock total (s) for the two canonical cases mirrored from
# scripts/profile-benchmarks.R: `scm` (run_scm) and `build_schedule`
# (run_scm refine_schedule = TRUE).
#
# Usage:
#   make compile                                   # build this branch
#   git worktree add -f /private/tmp/plant-other <ref> && ( cd /private/tmp/plant-other && make compile )
#   for r in 1 2; do
#     Rscript --no-init-file scripts/bench_ab.R "$(pwd)"            "this-r$r"
#     Rscript --no-init-file scripts/bench_ab.R /private/tmp/plant-other "other-r$r"
#   done
#
# Compare the RESULT lines. Trust same-session ratios, not absolute numbers
# across sessions (machine drift is large). See the profile-plant skill
# (.claude/skills/profile-plant/SKILL.md).

args  <- commandArgs(trailingOnly = TRUE)
path  <- if (length(args))     args[[1]] else "."
label <- if (length(args) > 1) args[[2]] else path
suppressMessages(pkgload::load_all(path, compile = FALSE, quiet = TRUE))

mk_scm <- function() expand_parameters(trait_matrix(0.0825, "lma"),
                                       scm_base_parameters("FF16"))
mk_bs  <- function() {
  p <- scm_base_parameters("FF16")
  p$strategies <- list(FF16_Strategy())
  p$birth_rate <- 0.1
  p
}
run_scm_case <- function() { run_scm(mk_scm()); invisible(NULL) }
run_bs_case  <- function() { run_scm(mk_bs(), refine_schedule = TRUE); invisible(NULL) }

med_ms    <- function(f, n = 7) {
  ts <- numeric(n)
  for (i in seq_len(n)) ts[i] <- system.time(f())[["elapsed"]]
  median(ts) * 1000
}
total20_s <- function(f) system.time(for (i in 1:20) f())[["elapsed"]]

invisible(run_scm_case()); invisible(run_bs_case())  # warm up (JIT/caches)

cat(sprintf("RESULT|%s|scm_one_ms=%.1f|bs_one_ms=%.1f|scm_20_s=%.3f|bs_20_s=%.3f\n",
            label, med_ms(run_scm_case), med_ms(run_bs_case),
            total20_s(run_scm_case), total20_s(run_bs_case)))
