# Same-machine A/B timing harness for the reverse-mode AD trait-gradient path
# (#472 scope B). Companion to scripts/bench_ab.R (SCM solve) and bench_tf24.R
# (TF24 hydraulics); this one times the GRADIENT calls.
#
# It times each gradient surface through its PUBLIC entry point only -- the public
# API now routes to the native (C++) entries, which harvest the frozen schedule
# from the live Patch and run the AD replay + reverse sweep in one fused call, so
# there is no separable R-side harvest left to time (the migration-era
# harvest_ms/impl_ms split it used to report is gone with the R harvest it
# measured). Per case it reports:
#   run_ms     -- the resident run_scm(save_RK45_cache=TRUE) (context; timed once)
#   public_ms  -- the end-to-end public gradient call (the headline)
#   val=<digest> -- a value digest of the Jacobian (same-session bit-identity)
#
# You cannot load two builds of `plant` in one R session, so build each and run
# this script once per build in a separate Rscript, interleaving to control for
# machine drift. Trust same-session ratios + the val digest (a build that changed
# any Jacobian moves it), not absolute numbers across sessions (see the
# profile-plant skill, .claude/skills/profile-plant/SKILL.md).
#
# Usage:
#   make compile                                   # build this branch
#   git worktree add -f /private/tmp/plant-other <ref> && ( cd /private/tmp/plant-other && make compile )
#   for r in 1 2; do
#     Rscript --no-init-file scripts/bench_gradient.R "$(pwd)"            "this-r$r"
#     Rscript --no-init-file scripts/bench_gradient.R /private/tmp/plant-other "other-r$r"
#   done
#
# Env vars:
#   BENCH_REPS   median over this many timed reps (default 7)
#   BENCH_CASES  comma-separated case subset (default all); e.g. "ff16_frozen,tf24f_census"

args  <- commandArgs(trailingOnly = TRUE)
path  <- if (length(args))     args[[1]] else "."
label <- if (length(args) > 1) args[[2]] else path
suppressMessages(pkgload::load_all(path, compile = FALSE, quiet = TRUE))

REPS <- as.integer(Sys.getenv("BENCH_REPS", "7"))

med_ms <- function(f, n = REPS) {
  ts <- numeric(n)
  for (i in seq_len(n)) ts[i] <- system.time(f())[["elapsed"]]
  median(ts) * 1000
}
# A cheap order-insensitive value digest: a Jacobian that changed anywhere moves it.
digest_jac <- function(x) {
  v <- if (is.list(x)) unlist(x[intersect(c("jacobian", "values", "d_state", "d_time",
                                            "time", "state", "gradient", "value"),
                                          names(x))], use.names = FALSE) else as.numeric(x)
  v <- v[is.finite(v)]
  sprintf("%.10g/%.10g/%d", sum(abs(v)), sum(v * seq_along(v)), length(v))
}

# ---- canonical cases (settings lifted verbatim from the test fixtures) --------

# FF16 frozen + offspring share one cached SCM (test-ff16-stand-gradient.R:7-12).
mk_ff16_basic <- function() {
  p <- scm_base_parameters("FF16")
  p <- add_strategies(p, trait_matrix(0.0825, "lma"), hyperpar = FF16_hyperpar,
                      birth_rate = list(20))
  p <- run_scm(p, Environment("FF16"), control(), refine_schedule = TRUE)$parameters
  run_scm(p, Environment("FF16"), control(save_RK45_cache = TRUE), refine_schedule = FALSE)
}
# FF16 single-species resident coupled (test-ff16-stand-gradient.R:66-73).
mk_ff16_resident <- function() {
  p <- scm_base_parameters("FF16")
  p <- add_strategies(p, trait_matrix(0.0825, "lma"), hyperpar = FF16_hyperpar,
                      birth_rate = list(20))
  p$node_schedule_times <- list(seq(0, 60, length.out = 11)); p$max_patch_lifetime <- 60
  run_scm(p, Environment("FF16"), control(save_RK45_cache = TRUE), refine_schedule = FALSE)
}
# FF16 two-species resident cross-species (test-ff16-stand-gradient.R:128-135).
mk_ff16_ms <- function() {
  p <- scm_base_parameters("FF16")
  p <- add_strategies(p, trait_matrix(c(0.0825, 0.2), "lma"), hyperpar = FF16_hyperpar,
                      birth_rate = list(20, 20))
  p$node_schedule_times <- list(seq(0, 70, length.out = 16), seq(0, 70, length.out = 16))
  p$max_patch_lifetime <- 70
  run_scm(p, Environment("FF16"), control(save_RK45_cache = TRUE), refine_schedule = FALSE)
}
# TF24f census + resident share one cached SCM (test-tf24f-census-gradient.R:20-28).
mk_tf24f <- function(H = 4L, n = 9L) {
  p <- scm_base_parameters("TF24f"); p$max_patch_lifetime <- H
  p <- add_strategies(p, trait_matrix(0.1978791, "lma"), hyperpar = TF24f_hyperpar,
                      birth_rate = list(20))
  p$node_schedule_times <- list(seq(0, H, length.out = n))
  ctlc <- control(shading_model = "crown-centre", GSS_tol_abs = 1e-9,
                  ode_tol_rel = 1e-4, ode_tol_abs = 1e-4, save_RK45_cache = TRUE)
  run_scm(p, Environment("TF24f"), ctlc, refine_schedule = FALSE)
}

# Each case: $run() builds the cached SCM / individual (timed once -> run_ms);
# $public(ctx, mets, tr) is the end-to-end public gradient call (its return value is
# digested into val). All public entries route to the native C++ path.
CASES <- list(
  ff16_frozen = list(
    run = mk_ff16_basic,
    mets = c("offspring_production", "LAI", "biomass", "size_moment"),
    tr = ff16_default_traits(),
    public = function(ctx, mets, tr) stand_gradient(ctx, metrics = mets, traits = tr)),

  ff16_offspring = list(
    run = mk_ff16_basic,
    mets = "offspring_production",
    tr = ff16_default_traits(),
    public = function(ctx, mets, tr) offspring_production_gradient(ctx, traits = tr)),

  ff16_resident_coupled = list(
    run = mk_ff16_resident,
    mets = c("LAI", "biomass", "size_moment"),
    tr = c("a_p1", "lma", "a_l1"),
    public = function(ctx, mets, tr) stand_gradient(ctx, metrics = mets, traits = tr,
                                                    feedback = "resident")),

  ff16_resident_ms = list(
    run = mk_ff16_ms,
    mets = c("LAI", "size_moment"),
    tr = c("lma", "a_p1"),
    public = function(ctx, mets, tr) stand_gradient(ctx, metrics = mets, traits = tr,
                                                    species = 1L, feedback = "resident")),

  tf24f_census = list(
    run = mk_tf24f,
    mets = c("LAI", "biomass", "size_moment"),
    tr = c("vcmax_25", "lma", "a_l1", "K_s"),    # TF24f pars subset (test-tf24f-census-gradient.R:75)
    public = function(ctx, mets, tr) tf24f_census_gradient_ad(ctx, metrics = mets, traits = tr)),

  tf24f_resident = list(
    run = mk_tf24f,
    mets = c("LAI", "size_moment"),
    tr = c("vcmax_25", "lma", "a_l1", "K_s"),    # TF24f pars subset (test-tf24f-census-gradient.R:132)
    public = function(ctx, mets, tr) tf24f_resident_census_gradient_ad(ctx, metrics = mets,
                                                                       traits = tr)),

  grow_individual = list(
    run = function() {                     # the "run" is the individual + env, not an SCM
      list(indv = Individual("FF16", "FF16_Env")(FF16_Strategy()), env = Environment("FF16"),
           targets = c(2, 5, 10))
    },
    mets = NA, tr = ff16_default_traits(),
    public = function(ctx, mets, tr) grow_individual_to_size_gradient(
      ctx$indv, ctx$targets, "height", ctx$env, traits = tr, time_max = 200))
)

# ---- run -----------------------------------------------------------------------

want <- Sys.getenv("BENCH_CASES", "")
sel  <- if (nzchar(want)) strsplit(want, ",")[[1]] else names(CASES)

for (nm in sel) {
  cs <- CASES[[nm]]
  if (is.null(cs)) { message("unknown case: ", nm); next }
  mets <- cs$mets; tr <- cs$tr

  run_ms <- med_ms(function() cs$run(), n = max(1L, REPS %/% 2L))   # solve is slow; fewer reps
  ctx <- cs$run()
  invisible(cs$public(ctx, mets, tr))                              # warm up
  public_ms <- med_ms(function() cs$public(ctx, mets, tr))
  val <- digest_jac(cs$public(ctx, mets, tr))

  cat(sprintf("RESULT|%s|%s|run_ms=%.3f|public_ms=%.3f|val=%s\n",
              label, nm, run_ms, public_ms, val))
}
