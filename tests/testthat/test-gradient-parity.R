# The parity gate: for every state the forward model returns a number for, the
# reverse returns a row or names a violated constraint.
#
# ⚠️ READ AS WRITTEN THIS IS SATISFIED BY A SWEEP THAT REFUSES EVERYTHING, which
# is why it is two checks and not one.
#
# The GATE is the first: nothing escapes unnamed. A raw error, a status outside
# the declared set, a finite-looking row under a refusal, an undefined metric
# that reads as a zero one -- each is a failure however few drivers reach it.
#
# The COVERAGE is the second, and it is a measurement rather than a threshold:
# which drivers answer, and for those that do not, the branch named. A regime
# that stops answering fails here; a regime that has never answered is listed
# below by name, so the difference between the two is a line of code rather than
# a reading of the output.

# The branches that have never returned a row, by the name the refusal carries.
# A refusal naming anything else is a regime that USED to answer, which is what
# this file exists to catch.
#
# ⚠️ ONE ENTRY, AND IT IS THE SWEEP'S OWN. `shade-death` left this list when the
# shut branches answered; the LIGHT FLOOR left it when the row it severs was
# found to be exactly zero for the model as evaluated -- every light below the
# floor gives a bit-identical census, so the honest row is the zero rather than a
# refusal, and what makes the zero readable is the clamp counter rather than a
# name here.
#
# What remains is not a branch of the model at all: the descent leaves the range
# a double holds. Nothing it computes is wrong -- see `ladder_range_refusal` --
# and it is a gap because the answer is unrepresentable rather than unknown.
parity_known_gaps <- ladder_range_refusal

# The drivers that reach it, by name. A count could not say this: a regime that
# stops answering and one that never answered both come back refused, and only
# the name separates them.
#
# Empty since the storage pool's charge and drain form came back: `shaded` and
# `clamped` sat here while the pool integrated past its own ceiling, and both
# answer now. Kept rather than deleted, because what it asserts now is that every
# driver answers, which is the reading a later refusal has to fail against.
parity_range_gap <- character(0)

# One driver, reduced to what the two directions each said.
#
# The trajectory is read AFTER the sweep and costs nothing: stand_gradient's own
# store_trajectory() has already kept it, where reading it first would run the
# model twice. Counted rather than carried -- one driver's states are 67 MB.
parity_of <- function(scm) {
  g <- stand_gradient(scm)
  refused_metric <- stand_gradient_refused(g)
  counts <- census_operating_point_counts_tf24(scm)[[1]]
  names(counts) <- census_operating_point_names_tf24()
  rec <- scm$store_trajectory()
  unrowed <- vapply(rec, function(r) sum(!is.finite(r$state)), 0)
  list(refused_metric = refused_metric, gradient = g,
       refused = any(refused_metric),
       reason = if (is.null(g$refusal[[1]])) NA_character_ else g$refusal[[1]]$reason,
       refusal = g$refusal[[1]],
       # Read after the sweep, so it is the sweep's own severances and not the
       # forward run's -- the two differ, and the forward one cannot say whether
       # a gradient carries a declared zero.
       clamp = census_clamp_counts_differentiated_tf24(scm)[[1]],
       kinds = counts[counts > 0],
       records = length(rec),
       unrowed_records = sum(unrowed > 0),
       unrowed_worst = max(unrowed))
}

# The regimes the reference capture uses, at the lifetime this file runs them
# for. Named once, in the helper: a driver listed here that the capture does not
# carry is a regime this file would report on and nothing would referee.
#
# `introductions` thins the default schedule, and it is a COST setting rather
# than a regime one -- the default derives its introduction count from the patch
# lifetime, so a driver pays for cohorts none of the checks below asks about.
# Measured over all five, run and swept, against the three axes this file
# reports: the verdict, the operating-point kinds reached, and the clamp sites
# the SWEEP met.
#
#   driver      steps 88 -> 20   seconds 88 -> 20   what moves at 20
#   wet            422 ->  368     16.4 ->   3.8    nothing
#   drought       2415 -> 1441     79.9 ->  12.1    nothing
#   shaded       11810 -> 3250    318.6 ->  21.5    nothing
#   clamped      11347 -> 3267    259.0 ->  17.6    nothing
#   seasonal      6582 -> 3723    156.1 ->  13.0    LOSES determined and
#                                                   hydraulic-shutdown
#
# ⚠️ SEASONAL KEEPS ITS WHOLE SCHEDULE AND IS THE REASON THIS IS PER DRIVER.
# It is this file's only route to a hydraulic shutdown, which is one of the
# operating-point kinds the gate exists to reach, and the thinned run reaches
# boundary-root-crit instead. Every other driver's verdict, kinds and clamp sites
# are the same at either length.
parity_introductions <- c(wet = 20L, drought = 20L, seasonal = NA_integer_,
                          shaded = 20L, clamped = 20L)

parity_drivers <- lapply(ladder_reference_regimes(), function(d) {
  keep <- parity_introductions[[d$name]]
  c(d, list(lifetime = 5, introductions = if (is.na(keep)) NULL else keep))
})

# Each driver is run and swept ONCE. The sweep is the whole cost here -- the run
# is free -- and three checks reading one sweep is the difference between a file
# that costs a minute and one that costs four.
parity_cache <- new.env(parent = emptyenv())

# What the answer is a function of, and nothing else: the compiled object, and the
# two files that define the fixture and reduce it. A key covering less than this
# could not see a rebuild, and a gate that cannot see a rebuild can report green
# for code that is gone -- which is the one failure a cache here must not have.
parity_key <- function() {
  dll <- getLoadedDLLs()[["plant"]][["path"]]
  paste(tools::md5sum(c(dll,
                        testthat::test_path("helper-gradient-ladder.R"),
                        testthat::test_path("test-gradient-parity.R"))),
        collapse = "-")
}

# Each regime is one forward run and one reverse sweep, and the five are
# independent, so they go out to separate processes: the file costs the slowest
# regime rather than their sum. Serial where forking is unavailable.
parity_compute <- function() {
  one <- function(d) {
    scm <- ladder_driver_stand(
      d$rain, d$lifetime,
      k_I = if (is.null(d$k_I)) 0.5 else d$k_I,
      amplitude = if (is.null(d$amplitude)) 0 else d$amplitude,
      introductions = d$introductions)
    c(list(name = d$name, census = stand_census(scm)), parity_of(scm))
  }
  n <- plant_test_cores(length(parity_drivers))
  if (.Platform$OS.type == "unix" && n > 1L) {
    parallel::mclapply(parity_drivers, one, mc.cores = n)
  } else {
    lapply(parity_drivers, one)
  }
}

parity_shared <- function() {
  if (!is.null(parity_cache$all)) {
    return(parity_cache$all)
  }
  # On disk only where a directory is named, so a plain run writes nothing and
  # the default behaviour is unchanged. The tiers run this file twice, and one
  # build's answer serves both.
  dir <- Sys.getenv("PLANT_TEST_CACHE", unset = "")
  file <- if (nzchar(dir)) file.path(dir, paste0("parity-", parity_key(), ".rds"))
          else ""
  if (nzchar(file) && file.exists(file)) {
    parity_cache$all <- readRDS(file)
    return(parity_cache$all)
  }
  parity_cache$all <- parity_compute()
  if (nzchar(file)) {
    dir.create(dir, showWarnings = FALSE, recursive = TRUE)
    saveRDS(parity_cache$all, file)
  }
  parity_cache$all
}

test_that("no driver the forward model answers leaves the reverse unnamed", {
  # The gate. Every clause fails on one driver, so none of them is a statement
  # about how many drivers there are.
  for (r in parity_shared()) {
    # The forward model has to have answered, or this driver says nothing about
    # the reverse and is not evidence either way.
    expect_true(all(is.finite(r$census)))

    entries <- r$gradient$gradient
    if (r$refused) {
      # Named, always. A refusal a reader cannot name is the raw error the
      # status channel replaced.
      expect_true(is.character(r$reason) && nchar(r$reason) > 0)

      # Located as far as the grain allows. One recording spans six stages and
      # every cohort in them, so a row it could not supply has nothing below the
      # species to name. The species survives only where the frame that raised it
      # knew one -- a refusal thrown from inside the recording does not.
      expect_true(r$refusal$species == -1 || r$refusal$species >= 1)
      if (r$refusal$species == -1) {
        message("      unlocated: raised inside the recording, which names no species")
      }
      # And an undefined metric must not read as a zero one -- the distinction
      # refusal exists for, asserted at the boundary rather than inside it. The
      # whole row goes, which is what makes one reason enough for it.
      expect_true(all(is.na(entries[r$refused_metric, , drop = FALSE])))
    } else {
      expect_true(all(is.finite(entries)))
      expect_gt(length(entries), 0)
    }
  }
})

test_that("every refusal names a branch that has never answered", {
  # The coverage half, and the reason it is not a count. A regime that stops
  # answering and one that never answered both come back "refused"; only the name
  # tells them apart, so the name is what is asserted.
  answered <- character(0)
  refused <- character(0)
  for (r in parity_shared()) {
    message(sprintf("  %-9s %-9s %s", r$name,
                    if (r$refused) "refused" else "answered",
                    paste(sprintf("%s=%.0f", names(r$kinds), r$kinds), collapse = " ")))
    if (r$refused) {
      refused <- c(refused, r$name)
      message(sprintf("      %s", substr(r$reason, 1, 140)))
      expect_true(any(vapply(parity_known_gaps,
                             function(g) grepl(g, r$reason, fixed = TRUE),
                             logical(1))))
    } else {
      answered <- c(answered, r$name)
    }
  }
  message(sprintf("  parity: %d of %d drivers answer; %s refuse",
                  length(answered), length(parity_shared()),
                  if (length(refused)) paste(refused, collapse = ", ") else "none"))

  # Asserted BY NAME in both directions, which is what keeps this from being a
  # count that a sweep answering by doing nothing would also satisfy: every
  # driver answers, and any that stops is named where it fails.
  expect_setequal(refused, parity_range_gap)
  expect_setequal(answered,
                  setdiff(vapply(parity_drivers, `[[`, "", "name"),
                          parity_range_gap))
})

test_that("the driver that reaches a clamp says so, answered or not", {
  # This file's non-vacuity guard, and it does not depend on the driver
  # answering: an answered gradient carrying a declared zero and one carrying no
  # clamp at all are the same numbers, so only the count separates them. Written
  # to hold whether or not a driver answers, because what it reads is the count a
  # descent left behind rather than the rows it returned.
  by_name <- stats::setNames(parity_shared(), vapply(parity_shared(),
                                                     function(r) r$name, ""))
  nm <- census_clamp_names_tf24()
  at <- function(driver, site) by_name[[driver]]$clamp[[match(site, nm)]]
  for (d in names(by_name)) {
    c_d <- by_name[[d]]$clamp
    message(sprintf("  %-9s %s", d,
                    paste(sprintf("%s=%.0f", nm[c_d > 0], c_d[c_d > 0]),
                          collapse = " ")))
  }

  # Every driver, including the control, severs its light-independent root rows:
  # above the rooting-depth cap the root profile stops reading height, and that
  # is the model rather than a guard. It is the non-vacuity guard because it
  # holds everywhere -- a build reporting zero here has lost the counter, not
  # found a cleaner stand.
  for (d in names(by_name)) {
    expect_gt(at(d, "rooting_depth"), 0)
  }

  # The light floor is the driver-specific one: the shaded driver reaches it and
  # the control does not, so these counts measure the driver rather than the clamp.
  # Both sites, because the crown one is where it binds first. Counted over the
  # steps the descent reached before it was refused, which is a subset of the run
  # and not none.
  expect_gt(at("clamped", "light_floor"), 0)
  expect_gt(at("clamped", "light_floor_crown"), 0)
  expect_gt(at("clamped", "light_floor_crown"), at("clamped", "light_floor"))
  expect_equal(at("wet", "light_floor"), 0)
  expect_equal(at("wet", "light_floor_crown"), 0)

  # And three sites fire on no driver at all, which is the guard census's own
  # entry and a different statement from "it held".
  for (s in c("soil_positivity", "rainfall", "infiltration")) {
    for (d in names(by_name)) {
      expect_equal(at(d, s), 0)
    }
  }
})

test_that("each driver reaches the branch it is here for", {
  # A fixture must be shown to reach what it tests. The drought and seasonal
  # drivers exist to put PINNED operating points on the gradient's path, and a
  # run that answered without ever reaching one would pass the gate while
  # testing nothing; the shaded driver exists to reach the one branch that has
  # no rows.
  by_name <- stats::setNames(parity_shared(), vapply(parity_shared(),
                                                     function(r) r$name, ""))
  reach <- function(nm, kind) {
    k <- by_name[[nm]]$kinds
    if (is.na(k[kind])) 0 else k[[kind]]
  }
  expect_gt(reach("drought", "boundary-crit"), 0)
  expect_gt(reach("seasonal", "boundary-crit"), 0)
  expect_gt(reach("shaded", "shade-death"), 0)
  expect_gt(reach("clamped", "shade-death"), 0)
  # And the control has to be a control: the wet driver never leaves the branch
  # the gradient was first built for, which is what makes it the one fixture a
  # regression shows up against cleanly.
  expect_equal(unname(reach("wet", "interior")), unname(sum(by_name[["wet"]]$kinds)))
})

test_that("no state the sweep replays is one a row cannot attach to", {
  # ⚠️ AN INFINITE STATE ENTRY IS NOT A FORWARD PROBLEM AND IS ALWAYS A REVERSE
  # ONE. Every reader of the cumulative hazard either tests is_finite or reads
  # exp(-mortality), so a run that parks one at +Inf completes with plausible
  # numbers and records a trajectory the sweep then replays; the gradient comes
  # back not-a-number for every metric with no refusal declared, which the gate
  # above cannot see because there is nothing to name. The introduction condition
  # was the source: -log(pr_estab) at an establishment probability of exactly
  # zero, which the shaded and clamped drivers reach and the other three do not.
  #
  # Asserted over the recorded states rather than over the gradient, because the
  # gradient is one number per column and says nothing about which state carried
  # it.
  for (r in parity_shared()) {
    message(sprintf("  %-9s %d records, %d carrying a state no row attaches to",
                    r$name, r$records, r$unrowed_records))
    expect_equal(r$unrowed_records, 0,
                 label = paste0(r$name, ": ", r$unrowed_records, " of ",
                                r$records, " recorded states carry up to ",
                                r$unrowed_worst, " non-finite entries"))
  }
})
