# Any strategy, in any environment, should run and return low growth and low
# fitness rather than an error. This file asks that of a Latin hypercube over
# ten traits and four environment axes, and holds every point to the four
# criteria plant-meta's model-robustness guidance names:
#
#   1. the run finishes;
#   2. growth and fitness are finite and non-negative;
#   3. no state leaves its domain;
#   4. fitness is monotone along named axes, on paired points.
#
# OPT-IN. Set PLANT_RUN_SWEEP=1 (and NOT_CRAN=true) to run it: it is meant for a
# nightly or pre-release run, not every push. It writes one row per point, with
# every axis's value, each criterion's verdict and every clamp site's count, to
# PLANT_SWEEP_OUT (a path; by default robustness-sweep.tsv in the session's
# temporary directory). The design is seeded, so a failing point is also
# `sweep_design()[point, ]`, and `sweep_point(model, that_row)` re-runs it alone.
#
# ⚠️ A FAILURE HERE IS A FINDING. Each one is either a model defect or an input
# outside the model's domain, and deciding which is the point of running it.
# Neither is answered by narrowing a range or relaxing a criterion until the
# file passes: a range is narrowed only with a reason that names the domain the
# model claims, written beside the range.
#
# WHAT EACH CRITERION READS.
#
# Growth is a cohort's height. TF24 and FF16 grow by allocating net production
# to new tissue, and net production may be negative (that is attainable, and the
# model represents it) but height may not fall. So every cohort's height is held
# non-decreasing across the recorded steps, and its height rate at the end is
# held finite and non-negative.
#
# The domain is checked on the states at every recorded step (each introduction
# and the end): every value finite; height positive; mortality, fecundity,
# heartwood, storage, survival-weighted offspring and density non-negative;
# patch survival in [0, 1]; soil moisture in [0, saturation]. The final ODE state
# must also pass the patch's own validity predicate. Between recorded steps the
# solver's own predicate rejects a step that leaves the non-negative states, so
# this is a check that it held, not the only guard.
#
# Monotonicity is asked only where it is physically certain, and so it is asked
# of mutants rather than of the environment. Moving rain, light or soil depth
# moves the resident's own shading and water use as well, and whether the
# resident's fitness then rises is an ecological question, not a property every
# correct model has. A mutant reads the field its resident recorded and does not
# move it, so a pair of mutants differing in one parameter differs in nothing
# else:
#
#   d_I, intrinsic mortality. It adds to the hazard at every age and enters no
#     growth rate, so a mutant with a higher d_I has the same fecundity schedule
#     and a strictly lower survival at every age past birth. Its fitness is
#     strictly lower wherever its twin's is positive.
#   a_d0, the establishment cost. It enters only the establishment probability,
#     1 / ((a_d0 * A_0 / P)^2 + 1) at the seed's net production P, which falls
#     as a_d0 rises. Its fitness is no higher; equality is allowed, because
#     where P is large the change rounds away.
#
# Clamp sites are counted on TF24's forward run (census_clamp_counts_tf24),
# every site the model and phylloptim's leaf list, phylloptim's vpd_leaf_floor
# among them. That site binds only with the energy balance on, which is why the
# third pass exists. FF16 has no counted clamp sites, so its rows carry none.

sweep_enabled <- function() identical(Sys.getenv("PLANT_RUN_SWEEP"), "1")

# The axes, each with its range and where the range comes from. `log = TRUE`
# samples uniformly in the logarithm, for an axis whose range spans more than a
# factor of about three. The curated sets cited are fieldtrial's
# (inst/sets/): alpine ash at Tumbarumba and mulga at Alice Springs, the two
# measured parameterisations the family has. A trait's range covers plant's
# default and both sets, with room either side.
sweep_axes <- function() {
  list(
    # Leaf mass per area [kg m^-2]. Default 0.198; alpine ash 0.206, mulga 0.458.
    lma      = list(lower = 0.05, upper = 0.6, log = TRUE),
    # Height at maturation [m]. Default 16.6; mulga 2.5, alpine ash 23.5.
    hmat     = list(lower = 1.5, upper = 40, log = TRUE),
    # Carboxylation capacity at 25 C [umol m^-2 s^-1]. Default 96; alpine ash
    # 81.5, mulga 174. jmax_25 moves with it at the default ratio, 1.64, so the
    # axis is photosynthetic capacity rather than an electron-transport
    # imbalance nobody measures.
    vcmax_25 = list(lower = 30, upper = 200, log = TRUE),
    # Stem water potential at 50% loss of conductivity [MPa, positive]. Default
    # 2.89; alpine ash 5.31, mulga 6.7 (both proxies in their sets).
    stem_P50 = list(lower = 1, upper = 8, log = TRUE),
    # Terminal-segment sapwood conductivity. Default 0.336; neither set has a
    # value, so a factor of five either side of the default.
    K_s      = list(lower = 0.336 / 5, upper = 0.336 * 5, log = TRUE),
    # Fine-root mass per leaf area [kg m^-2]. Default 0.07; no set has a value,
    # so a factor of four either side.
    a_r1     = list(lower = 0.07 / 4, upper = 0.07 * 4, log = TRUE),
    # Canopy light extinction coefficient. Default 0.5; 0.3 to 0.8 is the span
    # usually quoted for closed canopies.
    k_I      = list(lower = 0.3, upper = 0.8, log = FALSE),
    # Intrinsic mortality [yr^-1]. Default 0.01; from a long-lived tree to a
    # short-lived shrub.
    d_I      = list(lower = 0.002, upper = 0.1, log = TRUE),
    # Seed mass [kg]. Default 3.8e-5; alpine ash 1.2e-6, mulga 1.3e-5. 0.1 mg
    # to 1 g.
    omega    = list(lower = 1e-7, upper = 1e-3, log = TRUE),
    # Sapwood area per leaf area, the Huber value. Default 2.14e-4; alpine ash
    # 2.59e-4, mulga 1.57e-3.
    theta    = list(lower = 5e-5, upper = 2e-3, log = TRUE),
    # Mean rainfall [m yr^-1]. Default 1; the gradient ladder's drought and wet
    # regimes are 0.1 and 2.
    rain     = list(lower = 0.1, upper = 3, log = TRUE),
    # Seasonal amplitude: rainfall is rain * (1 + a sin(2 pi t)), so 0 is
    # constant and 1 is a dry season with none. The ladder's seasonal regime is
    # 1.
    season   = list(lower = 0, upper = 1, log = FALSE),
    # Shortwave above the canopy [W m^-2]. Default 900, a clear-sky midday
    # figure; down to a third of it.
    SWdown   = list(lower = 300, upper = 1000, log = FALSE),
    # Soil depth [m], in five equal layers. Default 1.5; mulga's set is 1.0, to
    # a hardpan.
    depth    = list(lower = 0.5, upper = 4, log = TRUE)
  )
}

# What each model has of the axes. FF16 carries seven of the traits (it has no
# photosynthetic capacity, stem vulnerability or stem conductivity parameter;
# its photosynthesis is a fixed curve of light) and none of the environment
# axes: its environment is the canopy's own light field and nothing else, with
# no water and no driver above the canopy. Its points are the same rows, read
# for the axes it has.
sweep_traits <- function(model) {
  all <- c("lma", "hmat", "vcmax_25", "stem_P50", "K_s", "a_r1", "k_I", "d_I",
           "omega", "theta")
  if (model == "FF16") setdiff(all, c("vcmax_25", "stem_P50", "K_s")) else all
}

# A Latin hypercube: each axis cut into n equal strata, one point per stratum,
# the strata paired at random across axes. Seeded, and the caller's random
# stream is left as it was.
sweep_design <- function(n = 200L, seed = 20261119L) {
  axes <- sweep_axes()
  u <- withr::with_seed(seed, vapply(seq_along(axes), function(i) {
    (sample.int(n) - stats::runif(n)) / n
  }, numeric(n)))
  out <- vapply(seq_along(axes), function(i) {
    a <- axes[[i]]
    if (a$log) {
      exp(log(a$lower) + u[, i] * (log(a$upper) - log(a$lower)))
    } else {
      a$lower + u[, i] * (a$upper - a$lower)
    }
  }, numeric(n))
  colnames(out) <- names(axes)
  data.frame(point = seq_len(n), out)
}

# The strategy's parameters at one point. The traits are set as parameters,
# with no trade-off function: TF24's trade-offs derive d_I and stem_P50 from
# wood density and K_s and refuse them as traits, and an axis that moved three
# parameters would make a failure three times as hard to attribute.
sweep_parameters <- function(model, x, energy_balance = FALSE) {
  p <- testing_parameters(model)
  p$max_patch_lifetime <- 10
  traits <- unlist(x[sweep_traits(model)])
  if (model == "TF24") {
    traits <- c(traits, jmax_25 = 1.64 * traits[["vcmax_25"]],
                use_energy_balance = as.numeric(energy_balance))
  }
  p <- add_strategies(p, trait_matrix(unname(traits), names(traits)),
                      hyperpar = NULL, birth_rate = 1)
  full <- p$node_schedule_times[[1]]
  p$node_schedule_times <- list(full[round(seq(1, length(full), length.out = 12))])
  p
}

sweep_environment <- function(model, x) {
  env <- Environment(model)
  if (model != "TF24") {
    return(env)
  }
  env$set_soil_layer_widths(rep(x$depth / 5, 5))
  env$set_soil_water_state(rep(0.5 * env$soil_moist_sat, 5))
  env$extrinsic_drivers_set_constant("SWdown", x$SWdown)
  t <- seq(0, 10, length.out = 481)
  env$extrinsic_drivers_set_variable(
    "rainfall", t, x$rain * (1 + x$season * sin(2 * pi * t)))
  env
}

# The mutants a point's monotonicity is read from: the resident's own strategy,
# then the same with d_I raised by 0.02 per year, then with a_d0 doubled. Built
# on the resident's strategy, so every other parameter is the point's, and on
# its introduction times.
sweep_mutants <- function(p) {
  base <- unlist(p$strategies[[1]]$pars[c("d_I", "a_d0")])
  m <- rbind(base, base + c(0.02, 0), base * c(1, 2))
  q <- p
  q$strategy_default <- p$strategies[[1]]
  add_mutant(q, trait_matrix(unname(m), colnames(m)), hyperpar = NULL,
             birth_rate = rep(1, 3))
}

# Every state check, as one string naming what failed and its worst value; ""
# where nothing did.
sweep_domain <- function(r, model) {
  td <- tidy(r)
  sp <- td$species
  sp <- sp[!is.na(sp$height), , drop = FALSE]
  bad <- character(0)
  worst <- function(x) signif(x[which.max(abs(x))], 3)
  states <- intersect(c("height", "mortality", "fecundity", "area_heartwood",
                        "mass_heartwood", "storage",
                        "offspring_produced_survival_weighted", "density"),
                      names(sp))
  for (s in states) {
    v <- sp[[s]]
    if (any(!is.finite(v))) {
      bad <- c(bad, sprintf("%s non-finite (%d)", s, sum(!is.finite(v))))
    } else if (s == "height" && any(v <= 0)) {
      bad <- c(bad, sprintf("height <= 0 (%s)", worst(v[v <= 0])))
    } else if (any(v < 0)) {
      bad <- c(bad, sprintf("%s < 0 (%s)", s, worst(v[v < 0])))
    }
  }
  pd <- td$steps$patch_density
  if (any(!is.finite(pd) | pd < 0 | pd > 1)) {
    bad <- c(bad, "patch survival outside [0, 1]")
  }
  if (model == "TF24") {
    theta <- td$env$soil_moist$soil_moist
    sat <- r$scenario$env$soil_moist_sat
    if (any(!is.finite(theta))) {
      bad <- c(bad, "soil moisture non-finite")
    } else {
      if (any(theta < 0)) bad <- c(bad, sprintf("soil moisture < 0 (%s)", worst(theta[theta < 0])))
      if (any(theta > sat)) {
        bad <- c(bad, sprintf("soil moisture above saturation by %s",
                              signif(max(theta - sat), 3)))
      }
    }
  }
  patch <- r$solver$patch
  y <- patch$ode_state
  if (any(!is.finite(y)) || !isTRUE(patch$ode_state_valid(y))) {
    bad <- c(bad, "final ODE state refused by the patch's validity predicate")
  }
  paste(bad, collapse = "; ")
}

# Growth: every cohort's height non-decreasing over the recorded steps, and its
# height rate at the end finite and non-negative. "" where both hold.
sweep_growth <- function(r) {
  sp <- tidy(r)$species
  sp <- sp[!is.na(sp$height), , drop = FALSE]
  bad <- character(0)
  drop <- tapply(sp$height, sp$node, function(h) {
    if (length(h) < 2L) 0 else min(diff(h))
  })
  if (any(drop < 0)) {
    bad <- c(bad, sprintf("height fell (by %s)", signif(-min(drop), 3)))
  }
  rates <- vapply(r$solver$patch$species[[1]]$nodes,
                  function(nd) nd$ode_rates[[1]], 0)
  if (any(!is.finite(rates))) {
    bad <- c(bad, "height rate non-finite")
  } else if (any(rates < 0)) {
    bad <- c(bad, sprintf("height rate < 0 (%s)", signif(min(rates), 3)))
  }
  paste(bad, collapse = "; ")
}

sweep_clamp_names <- function() census_clamp_names_tf24()

# One point, run and judged. Never throws: a point that errors is a row saying
# so, so the sweep reports every point rather than the first to fail.
sweep_point <- function(model, x, energy_balance = FALSE) {
  started <- proc.time()[["elapsed"]]
  row <- list(finished = FALSE, error = "", fitness = NA_real_, growth = "",
              domain = "", fitness_d_I = NA_real_, fitness_a_d0 = NA_real_,
              fitness_twin = NA_real_, monotone = "")
  clamps <- stats::setNames(rep(NA_real_, length(sweep_clamp_names())),
                            paste0("clamp_", sweep_clamp_names()))
  got <- tryCatch({
    p <- sweep_parameters(model, x, energy_balance)
    r <- run(p, sweep_environment(model, x), record = "steps")
    row$finished <- isTRUE(r$complete)
    row$fitness <- r$net_reproduction_ratios[[1]]
    if (model == "TF24") {
      clamps[] <- census_clamp_counts_tf24(r$solver)[[1]]
    }
    row$growth <- sweep_growth(r)
    row$domain <- sweep_domain(r, model)
    # Last, because it moves the resident's solver on.
    f <- run_mutant(r, sweep_mutants(p))$net_reproduction_ratios
    row$fitness_twin <- f[[1]]
    row$fitness_d_I <- f[[2]]
    row$fitness_a_d0 <- f[[3]]
    mono <- character(0)
    if (!all(is.finite(f))) {
      mono <- c(mono, "mutant fitness non-finite")
    } else {
      if (if (f[[1]] > 0) f[[2]] >= f[[1]] else f[[2]] != 0) {
        mono <- c(mono, sprintf("d_I + 0.02 did not lower fitness (%s -> %s)",
                                signif(f[[1]], 4), signif(f[[2]], 4)))
      }
      if (f[[3]] > f[[1]]) {
        mono <- c(mono, sprintf("2 a_d0 raised fitness (%s -> %s)",
                                signif(f[[1]], 4), signif(f[[3]], 4)))
      }
    }
    row$monotone <- paste(mono, collapse = "; ")
    NULL
  }, error = function(e) conditionMessage(e))
  if (!is.null(got)) {
    row$error <- gsub("[\t\n]+", " ", got)
  }
  data.frame(row, as.list(clamps), seconds = proc.time()[["elapsed"]] - started,
             stringsAsFactors = FALSE, check.names = FALSE)
}

# Which of the four criteria a row fails, by name.
sweep_failures <- function(res) {
  fit_ok <- is.finite(res$fitness) & res$fitness >= 0
  data.frame(
    finishes = !res$finished | nzchar(res$error),
    finite_non_negative = res$finished & !nzchar(res$error) &
      (!fit_ok | nzchar(res$growth)),
    in_domain = nzchar(res$domain),
    monotone = nzchar(res$monotone))
}

sweep_run <- function(model, design, energy_balance = FALSE) {
  rows <- lapply(seq_len(nrow(design)), function(i) design[i, , drop = FALSE])
  one <- function(x) {
    cbind(x, sweep_point(model, x, energy_balance))
  }
  n <- plant_test_cores(length(rows))
  out <- if (.Platform$OS.type == "unix" && n > 1L) {
    parallel::mclapply(rows, one, mc.cores = n, mc.preschedule = FALSE)
  } else {
    lapply(rows, one)
  }
  # A worker that died returns a condition rather than a row; it is reported as
  # the point not finishing, with the worker's message, rather than lost.
  out <- lapply(seq_along(out), function(i) {
    if (inherits(out[[i]], "try-error") || inherits(out[[i]], "condition")) {
      cbind(rows[[i]], sweep_point_failed(model, as.character(out[[i]])))
    } else {
      out[[i]]
    }
  })
  res <- do.call(rbind, out)
  res$set <- if (energy_balance) "TF24 energy balance" else model
  res
}

sweep_point_failed <- function(model, message) {
  clamps <- stats::setNames(rep(NA_real_, length(sweep_clamp_names())),
                            paste0("clamp_", sweep_clamp_names()))
  data.frame(finished = FALSE, error = gsub("[\t\n]+", " ", message),
             fitness = NA_real_, growth = "", domain = "",
             fitness_d_I = NA_real_, fitness_a_d0 = NA_real_,
             fitness_twin = NA_real_, monotone = "", as.list(clamps),
             seconds = NA_real_, stringsAsFactors = FALSE, check.names = FALSE)
}

# The failing rows as a table a reader can act on: the point, the criteria it
# failed, and what each said.
sweep_report <- function(res) {
  fail <- sweep_failures(res)
  hit <- rowSums(fail) > 0
  if (!any(hit)) {
    return(invisible(NULL))
  }
  what <- apply(fail[hit, , drop = FALSE], 1, function(z) {
    paste(names(z)[z], collapse = ",")
  })
  detail <- with(res[hit, , drop = FALSE],
                 trimws(paste(error, growth, domain, monotone)))
  tab <- data.frame(set = res$set[hit], point = res$point[hit],
                    failed = what, detail = substr(detail, 1, 160))
  message(paste(utils::capture.output(print(tab, row.names = FALSE)),
                collapse = "\n"))
}

test_that("every point of a wide trait and environment sweep meets the four criteria", {
  skip_if_not(sweep_enabled(),
              "the robustness sweep is opt-in: set PLANT_RUN_SWEEP=1 (nightly or before a release)")
  skip_on_cran()

  design <- sweep_design()
  started <- proc.time()[["elapsed"]]
  res <- rbind(sweep_run("FF16", design),
               sweep_run("TF24", design),
               sweep_run("TF24", design, energy_balance = TRUE))
  elapsed <- proc.time()[["elapsed"]] - started

  out <- Sys.getenv("PLANT_SWEEP_OUT",
                    unset = file.path(tempdir(), "robustness-sweep.tsv"))
  utils::write.table(res, out, sep = "\t", quote = FALSE, row.names = FALSE)
  message(sprintf("robustness sweep: %d points in %.0f s, written to %s",
                  nrow(res), elapsed, out))
  sweep_report(res)

  # Every clamp site's incidence, summed over the points of each TF24 pass, so
  # the run says where each binds outside the ladder's stands.
  clamp_columns <- paste0("clamp_", sweep_clamp_names())
  # Counted over the points that finished: a point that errored has no count.
  for (set in unique(res$set[res$set != "FF16"])) {
    counts <- res[res$set == set & !is.na(res$clamp_rooting_depth),
                  clamp_columns, drop = FALSE]
    message(sprintf("  %-20s %s", set, paste(sprintf(
      "%s %d/%d", sweep_clamp_names(), colSums(counts > 0), nrow(counts)),
      collapse = ", ")))
  }

  fail <- sweep_failures(res)
  for (set in unique(res$set)) {
    mine <- fail[res$set == set, , drop = FALSE]
    for (criterion in names(fail)) {
      expect_equal(sum(mine[[criterion]]), 0L,
                   label = sprintf("%s: points failing '%s' (of %d)", set,
                                   criterion, nrow(mine)))
    }
  }
})
