
source(
  file.path(
    rprojroot::find_testthat_root_file(),
    "FF16_reference", "make_reference_plant.R")
)

# Edition-3 replacement for the deprecated `expect_is()`. Preserves the exact
# `inherits()` semantics expect_is relied on (works for S3/R6/RcppR6 classes as
# well as base types such as "list"/"numeric"/"matrix"), while giving a useful
# failure message. Using this avoids the per-call-site ambiguity of choosing
# between expect_s3_class / expect_type / expect_s4_class across ~70 sites.
expect_inherits <- function(object, class, info = NULL, label = NULL) {
  act <- testthat::quasi_label(rlang::enquo(object), label, arg = "object")
  testthat::expect(
    inherits(act$val, class),
    sprintf("%s does not inherit from `%s`; its class is `%s`.",
            act$lab, paste(class, collapse = "/"),
            paste(class(act$val), collapse = "/")),
    info = info
  )
  invisible(act$val)
}

# Compare two `Internals` reference objects by value. Edition-3 waldo recurses
# into the underlying C++ `.ptr`, which differs between two equivalent-but-distinct
# objects; comparing the exposed numeric bindings instead asserts equality of the
# *state* (which is what edition 2 effectively compared for these reference types).
expect_equal_internals <- function(object, expected, ...) {
  to_state <- function(x) list(states = x$states, rates = x$rates, auxs = x$auxs)
  expect_equal(to_state(object), to_state(expected), ...)
}

# One compute_rates call on a standalone TF24 individual, returning its named
# aux vector. The cheapest probe of the hydraulic path there is: no SCM, no
# patch, one leaf optimisation. Shared by test-tf24-root-pars.R and
# test-tf24-stem-hydraulics.R.
#
# `soil_theta` is soil water content per layer (m^3 m^-3) -- note that this is
# an entirely different `theta` from the strategy's `pars$theta`, which is leaf
# area per sapwood area. Both are reachable in the same session; the argument is
# named for the one it is.
tf24_probe <- function(strategy, soil_theta, height = 5) {
  ind <- TF24_Individual(strategy)
  ind$set_state("height", height)
  env <- TF24_Environment()
  env$set_soil_number_of_depths(length(soil_theta))
  env$set_soil_water_state(soil_theta)
  ind$compute_rates(env)
  stats::setNames(ind$internals$auxs, ind$aux_names)
}

# plant's testing parameter set, quietly: today's C++ defaults with the frozen
# trade-off functions attached (see ?default_parameters). Every test that builds
# strategies from traits uses this rather than scm_base_parameters(), which
# attaches no trade-offs, so the goldens and baselines keep their numbers.
testing_parameters <- function(model) {
  default_parameters(model, quiet = TRUE)
}

# Every model plant ships, for tests that check something of all of them.
all_models <- c("FF16", "K93", "TF24", "TF24f", "TF24t", "TF24w")

# One FF16 species at the testing parameters, run for `lifetime` years on the
# default node schedule.
ff16_parameters <- function(lifetime) {
  p <- testing_parameters("FF16")
  p$max_patch_lifetime <- lifetime
  p$node_schedule_times <- list()
  add_strategies(p, trait_matrix(0.0825, "lma"))
}

# d(metric)/d(name) for species `species` by a centred difference of two runs at
# x * (1 +- rel), where `measure(p)` runs p and returns its named metrics.
central_difference <- function(p, species, name, metric, measure, rel = 1e-4) {
  x <- p$strategies[[species]]$pars[[name]]
  h <- rel * abs(x)
  at <- function(v) {
    q <- p
    q$strategies[[species]]$pars[[name]] <- v
    measure(q)[[metric]]
  }
  (at(x + h) - at(x - h)) / (2 * h)
}

# Soil layer widths grading geometrically from a thin surface layer, summing to
# `depth` to the bit. The generator moved to fieldtrial (soil_widths_graded());
# the tests that probe a graded column keep this copy of its arithmetic.
graded_widths <- function(depth = 1.5, n = 5, top = 0.02) {
  n <- as.integer(n)
  stopifnot(n >= 1L, top < depth / n)
  if (n == 1L) {
    return(depth)
  }
  total <- function(r) top * (r^n - 1) / (r - 1)
  upper <- 2
  while (total(upper) < depth) {
    upper <- upper * 2
  }
  r <- stats::uniroot(function(r) total(r) - depth,
                      c(1 + 1e-12, upper), tol = .Machine$double.eps^0.75)$root
  out <- top * r^(seq_len(n) - 1L)
  out[[n]] <- depth - sum(out[-n])
  out
}

# Rainfall as Poisson pulses over (from, to]: exponential inter-arrival at
# `rate` per year, exponential depths with mean `mean_amount`. The generator
# moved to fieldtrial (pulse_series_poisson()); TF24w's stand keeps its
# arithmetic so the seeded series the test was tuned on is unchanged.
poisson_pulses <- function(from, to, rate, mean_amount) {
  n <- stats::rpois(1, rate * (to - from))
  time <- sort(stats::runif(n, from, to))
  amount <- stats::rexp(n, rate = 1 / mean_amount)
  keep <- time > from
  data.frame(time = time[keep], amount = amount[keep])
}

# NOTE: TF24f (#525) is intentionally excluded from these generic cross-strategy
# lists. It is an experimental variant with a 6th ODE state (tracked root-collar
# psi) and relies on birth-initialisation via Node::compute_initial_conditions, so
# it does not satisfy the generic invariants these lists drive (5-state layout;
# growth of a standalone Individual stepped a few times, which bypasses that
# initialisation). TF24f is covered via the SCM in test-strategy-tf24f.R. TF24t
# is excluded for the same reasons (four appended states, seeded in
# set_initial_states) and is covered in test-strategy-tf24t.R.

get_list_of_strategy_types <- function() {
  list(
    FF16=FF16_Strategy,
    TF24=TF24_Strategy,
    K93=K93_Strategy
    )
}

get_list_of_environment_types <- function() {
  list(
    FF16="FF16_Env",
    TF24="TF24_Env",
    K93="K93_Env"
    )
}

# How many processes a test may fork.
#
# ⚠️ TWO UNDER `R CMD check`, AND THAT IS NOT A PREFERENCE. The check sets
# `_R_CHECK_LIMIT_CORES_`, and `mclapply` then ERRORS rather than warning on a
# third process -- "3 simultaneous processes spawned" -- so `detectCores() - 1`
# fails the file outright on any runner with more than three cores. plant's
# workflow runs `--as-cran` with `error-on: "error"`, so that is a red leg on
# every operating system.
#
# Named once because two files fork and a rule kept in both is a rule that has to
# stay in sync; the third file to fork will find it here.
plant_test_cores <- function(jobs) {
  if (nzchar(Sys.getenv("_R_CHECK_LIMIT_CORES_"))) {
    return(min(jobs, 2L))
  }
  min(jobs, max(1L, parallel::detectCores() - 1L))
}
