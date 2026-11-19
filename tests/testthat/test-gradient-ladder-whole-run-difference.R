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
# The pinned grid is the run's own, so the reference is the derivative of the
# DISCRETISED run, and any change that moves the steps a run takes ages it.
# Re-capturing then is legitimate only once a difference with every step pinned
# agrees with the sweep on the rows that moved: that is what separates an aged
# reference from a wrong gradient.
#
# A residual is normalised per metric, over the columns the reference carries for
# it, which is the convention `ladder_matrix_residual` uses. Normalising per
# column instead makes every column whose whole row is at the difference's own
# floor look infinitely wrong: `a_f2` reads 4e-15 against an exact zero, and a
# ratio says that is a total disagreement where the honest reading is that both
# are zero. So each row the reference resolves is also held to its own scale,
# which the per-metric residual alone cannot see.

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

# The columns this reference carries no row for, and the reason is not that
# nobody captured them. Two default to 0, so the capture's RELATIVE step
# `abs(value) * rel` is 0 and the difference moves no parameter at all; their
# rows record that refusal and carry no metric, which reference_rows() drops.
# S_D is captured: it enters no rate, so its census rows are zero by
# construction (ladder_zero_outside_every_rate), and the sweep carries it as a
# fitness column.
#
# ⚠️ ASSERTED BOTH WAYS: `reference_compare` drops a reference row whose column
# the sweep does not carry -- which is right, names change -- and has no way to
# notice a column the REFERENCE does not carry, so a reference missing columns
# would pass with nothing said.
reference_uncaptured_columns <- function() {
  c("TF24_floor_lambda_o", "recruitment_decay")
}

reference_compare <- function(regime, rows) {
  solver <- reference_stand(regime)
  got <- stand_gradient(solver, hyperpar = NULL)
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
       # Which branches the run reached, summed over species, for whoever reads a
       # failure.
       kinds = stats::setNames(Reduce(`+`, census_operating_point_counts_tf24(solver)),
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
       spread_raw = mine$spread,
       residual = abs(observed - mine$converged) /
         unname(scale[mine$metric]))
}

# The regimes whose sweep leaves the range a double holds, so this reference
# has nothing to referee on them. Named rather than inferred: a regime that stops
# answering and one that never answered both arrive refused.
#
# Empty, and kept rather than deleted, because what it asserts is that nothing
# is refused.
reference_range_gap <- character(0)

test_that("the sweep agrees with a difference of whole runs, over five regimes", {
  rows <- reference_rows()
  regimes <- ladder_reference_regimes()
  expect_setequal(unique(rows$regime), vapply(regimes, `[[`, "", "name"))

  # The reference covers every column the sweep carries, bar the declared one.
  # Asserted before any comparison, because a missing row is not a disagreement
  # -- it is a column nothing looked at, and the rung passes either way.
  carried <- setdiff(names(TF24_Model()$pars),
                     names(census_undifferentiable("TF24")))
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
    # reads 8.6e-04, on clamped's `theta`, whose four steps themselves spread by
    # as much; no other regime reads above 3e-07. `theta` and `omega` reach the
    # census through the channels the leaf boundary carries and are the last
    # columns to resolve.
    tolerance <- pmax(3 * r$spread, 2e-3)
    over <- live & r$residual > tolerance
    # ⚠️ THIS RUNG IS THE ONLY ONE THAT SEES THE KNOT GRID, through this floor.
    # The light field is tabulated against knots at u_k * height_max, and the
    # canopy top is a cohort's height, so the grid moves with it. With the knot
    # positions passivised the gradient carries only the values and the slopes,
    # and `shaded` and `clamped` read 9.9e-03 and 8.0e-03 here. Every other rung
    # takes the field AS values and slopes, so none of them can see it. Widening
    # this floor to pass would bury it.
    expect_equal(sum(over), 0L,
                 label = paste0(r$name, ": ", sum(over), " column(s) past the ",
                                "reference's own spread, worst ",
                                r$column[[which.max(ifelse(live, r$residual, 0))]],
                                " at ",
                                signif(max(r$residual[live]), 3)))

    # And each row the reference resolves, on its own scale. The check above
    # divides by the metric's largest entry, so a small row passes it whatever
    # the sweep says: a zero passes most of them. A row whose four steps agree to
    # 1e-4 of its value is held to 1e-4 of it, or to three times its spread. The
    # floor is the reference's own resolution: it reads a column at one step for
    # every metric, so a sharply curved column carries that step's truncation
    # error. The row nearest its tolerance, shaded's leaf_area in 2.a_f2, is read
    # at 1e-3 and sits 7e-05 from the sweep, which a Richardson extrapolation of
    # the difference matches to 1e-05: the margin there is the reference's error,
    # not the sweep's. A reference aged by a moved step grid fails this floor
    # first.
    resolved <- live & r$spread_raw < 1e-4 * abs(r$reference)
    expect_gt(sum(resolved), 150)
    own <- pmax(3 * r$spread_raw, 1e-4 * abs(r$reference))
    gap <- abs(r$observed - r$reference)
    off <- resolved & gap > own
    expect_equal(sum(off), 0L,
                 label = paste0(r$name, ": ", sum(off), " resolved row(s) past ",
                                "their own scale, worst ",
                                paste(r$metric, r$column)[[which.max(
                                  ifelse(resolved, gap / own, 0))]],
                                " at ", signif(max((gap / own)[resolved]), 3),
                                " of its tolerance"))
  }

})
