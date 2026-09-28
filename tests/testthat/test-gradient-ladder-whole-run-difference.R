# The captured reference, which is the one instrument that shares no arithmetic
# with the sweep.
#
# Every other rung checks the gradient against something the same code computes a
# second way: a forward tangent of the same recording, a rank structure, an
# identity between two assemblies. Those catch an assembly error and cannot catch
# a wrong row, because a wrong row is wrong in both directions. This one compares
# against a difference of whole runs at four step sizes, taken on a stand where
# two species compete, and it is what a change from a differenced row to a
# closed-form one has to answer to.
#
# ⚠️ THE REFERENCE IS DATA AND IS NEVER REWRITTEN BY A TEST. It lives in
# reference/reference-gradient.tsv and `scripts/capture-reference-gradient.R`
# writes it, which is a thing to run when the INSTRUMENT moves -- the step
# ladder, the pinned grid, the parameter set -- and never to make a run green. A
# reference regenerated from the model it referees is not a reference, so a
# disagreement here is a finding.
#
# A residual is normalised per metric, over the columns the reference carries for
# it, which is the convention `ladder_matrix_residual` uses. Normalising per
# column instead makes every column whose whole row is at the difference's own
# floor look infinitely wrong: `a_f2` reads 4e-15 against an exact zero, and a
# ratio says that is a total disagreement where the honest reading is that both
# are zero.

reference_rows <- function() {
  path <- testthat::test_path("reference", "reference-gradient.tsv")
  skip_if_not(file.exists(path), "no captured reference")
  # The numbers are written right-padded, so they arrive as character unless the
  # padding is stripped -- and a silently character column makes every comparison
  # below an error rather than a disagreement.
  rows <- utils::read.delim(path, stringsAsFactors = FALSE, strip.white = TRUE)
  for (column in c("converged", "spread", "step")) {
    rows[[column]] <- as.numeric(rows[[column]])
  }
  # A row with no metric is the capture's own record of a column it could not
  # difference, and carries no number to referee against.
  rows[!is.na(rows$metric) & nzchar(rows$metric) & !is.na(rows$converged), ,
       drop = FALSE]
}

# The stand the reference was captured on: two species that compete, at the
# introduction times the capture used, under one regime. Built from the same
# helpers the capture used, so a regime named in both places is one stand.
reference_stand <- function(regime) {
  p <- ladder_parameters(c("fast", "slow"), k_I = regime$k_I)
  p$node_schedule_times <- list(c(0, 0.63), c(0, 0.41))
  ladder_run(p, ladder_control(),
             env = ladder_environment(regime$rain,
                                      if (is.null(regime$amplitude)) 0
                                      else regime$amplitude))
}

# The two columns this reference carries no row for, and the reason is not that
# nobody captured them: both default to 0, so the capture's RELATIVE step
# `abs(value) * rel` is 0 and the difference moves no parameter at all. Their
# rows record that refusal and carry no metric, which reference_rows() drops.
#
# ⚠️ ASSERTED BOTH WAYS, because the alternative is what this rung did while
# passing at 13: `reference_compare` drops a reference row whose column the sweep
# does not carry -- which is right, names change -- and has no way to notice a
# column the REFERENCE does not carry. Forty of forty-eight were refereed and
# nothing said so.
reference_uncaptured_columns <- function() {
  c("TF24_floor_lambda_o", "recruitment_decay")
}

reference_compare <- function(regime, rows) {
  scm <- reference_stand(regime)
  got <- stand_gradient(scm)
  mine <- rows[rows$regime == regime$name, , drop = FALSE]
  column <- paste0(mine$species, ".", mine$parameter)
  # A column the sweep does not carry is the capture's business rather than this
  # check's: it recorded every registered parameter, and the sweep reports a
  # class for the ones it cannot answer.
  keep <- column %in% colnames(got$gradient) &
    mine$metric %in% rownames(got$gradient)
  mine <- mine[keep, , drop = FALSE]
  column <- column[keep]
  observed <- got$gradient[cbind(mine$metric, column)]
  scale <- tapply(pmax(abs(mine$converged), abs(observed)), mine$metric, max)
  list(name = regime$name, column = column, metric = mine$metric,
       parameter = mine$parameter,
       observed = observed,
       # Which branches the run reached, summed over species, because the
       # tolerance below depends on one of them rather than on the regime's name.
       kinds = stats::setNames(Reduce(`+`, census_operating_point_counts_tf24(scm)),
                               census_operating_point_names_tf24()),
       # Refusal is metric-level, so it broadcasts to the rows this picked out.
       refused = unname(stand_gradient_refused(got)[mine$metric]),
       # And the reason, so a regime refusing every metric can be checked
       # against the gap it is meant to be rather than only counted.
       reason = if (any(stand_gradient_refused(got)))
                  got$refusal[[which(stand_gradient_refused(got))[[1]]]]$reason
                else NA_character_,
       reference = mine$converged,
       # The capture's own resolution and this rung's disagreement, divided by
       # one scale in one place: the capture reports the gap between the two
       # readings it chose between in the metric's own units, and normalising it
       # anywhere but here would be a second division that hides this one.
       spread = mine$spread / unname(scale[mine$metric]),
       residual = abs(observed - mine$converged) /
         unname(scale[mine$metric]))
}

# The regimes whose descent leaves the range a double holds, so this reference
# has nothing to referee on them. Named rather than inferred: a regime that stops
# answering and one that never answered both arrive refused.
#
# Empty since the storage pool's charge and drain form came back: `shaded` and
# `clamped` sat here while the pool integrated past its own ceiling, and the
# reference has carried 280 rows for each of them all along. Kept rather than
# deleted, because what it asserts now is that nothing is refused.
reference_range_gap <- character(0)

test_that("the sweep agrees with a difference of whole runs, over five regimes", {
  rows <- reference_rows()
  regimes <- ladder_reference_regimes()
  expect_setequal(unique(rows$regime), vapply(regimes, `[[`, "", "name"))

  # The reference covers every column the sweep carries, bar the declared one.
  # Asserted before any comparison, because a missing row is not a disagreement
  # -- it is a column nothing looked at, and the rung passes either way.
  carried <- setdiff(names(TF24_Strategy()$pars),
                     names(census_undifferentiable_tf24()))
  expect_setequal(setdiff(carried, unique(rows$parameter)),
                  reference_uncaptured_columns())

  # A regime is one run and one sweep and the five are independent, so they go to
  # separate processes -- but a forked worker that fails reports only that all
  # cores encountered errors, so each carries its own condition back.
  n <- plant_test_cores(length(regimes))
  guarded <- function(regime) {
    tryCatch(reference_compare(regime, rows), condition = function(e) e)
  }
  results <- if (.Platform$OS.type == "unix" && n > 1L) {
    parallel::mclapply(regimes, guarded, mc.cores = n)
  } else {
    lapply(regimes, guarded)
  }
  for (r in results) {
    if (inherits(r, "condition")) {
      stop(conditionMessage(r), call. = FALSE)
    }
  }

  # A regime refusing every metric has no column to referee, and asserting the
  # residual over none of them would pass vacuously. Checked by NAME against the
  # gap it is meant to be, then set aside.
  for (r in results) {
    if (all(r$refused)) {
      message(sprintf("  %-9s refused: %s", r$name, substr(r$reason, 1, 120)))
      expect_true(r$name %in% reference_range_gap)
      expect_true(grepl(ladder_range_refusal, r$reason, fixed = TRUE))
    }
  }
  names_of <- vapply(results, function(r) r$name, "")
  expect_setequal(names_of[vapply(results, function(r) all(r$refused), TRUE)],
                  reference_range_gap)

  for (r in results) {
    if (all(r$refused)) {
      next
    }
    # Answered columns only. Whether a refused metric is refused for the right
    # reason is the parity and sweep rungs' question. The declared zeros are out
    # for a separate reason: the declared-zero file referees them against the
    # model's own claim, and a relative residual has nothing to say about a
    # column that is exactly zero against a difference that is merely small.
    live <- !r$refused & !(r$parameter %in% ladder_zero_by_construction())
    expect_gt(sum(live), 200)
    # The reference's own resolution, with a floor: where its four steps agreed
    # to round-off, the sweep is still only asked to agree to the truncation the
    # coarsest of them carries. Over 270 answered columns a regime, the worst
    # reads 1.1e-03 on drought, 7.6e-04 on seasonal and 6.3e-05 on wet -- and all
    # three are `theta` or `omega`, which reach the census through the channels
    # the leaf boundary carries and are the last columns to resolve.
    tolerance <- pmax(3 * r$spread, 2e-3)
    over <- live & r$residual > tolerance
    # ⚠️ THIS RUNG IS WHAT FOUND THE KNOT GRID, and the floor is what found it.
    # `shaded` and `clamped` read 9.9e-03 and 8.0e-03 here while the light field
    # was held against knots at u_k * height_max whose positions were passivised:
    # the canopy top is a cohort's height, so the grid moved with it and the
    # gradient carried only the values and the slopes. Every other rung takes the
    # field AS values and slopes, so none of them could see it. Widening this
    # floor to pass would have buried it.
    expect_equal(sum(over), 0L,
                 label = paste0(r$name, ": ", sum(over), " column(s) past the ",
                                "reference's own spread, worst ",
                                r$column[[which.max(ifelse(live, r$residual, 0))]],
                                " at ",
                                signif(max(r$residual[live]), 3)))
  }

})
