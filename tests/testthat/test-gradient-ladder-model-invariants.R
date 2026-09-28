# The floor: what holds before any reference exists.
#
# These need no reference at all and are gated behind nothing. They are the
# cheapest assurance available and they run at every rung after this one. If any
# of them fails, nothing above is worth building yet.
#
# Tolerance is exactly zero for the decomposition family. There is no margin, so
# their whole assurance is the non-vacuity clause: prove the two paths computed
# something before believing they agree.

# The objective itself.

test_that("the fixtures sit in the regime their checks are declared for", {
  patch <- ladder_patch_two_by_two()
  report <- ladder_regime_report(patch, "patch")
  message("\npatch fixture regime:")
  for (i in seq_len(nrow(report))) {
    message(sprintf("  %-58s %-5s %s", report$assertion[[i]],
                    report$ok[[i]], report$value[[i]]))
  }
  expect_true(all(report$ok))
  # Crossing is asserted by the check that needs it rather than by the table,
  # because the fixtures come both crossed and not and a failure has to be
  # localisable to the crossing.
  expect_equal(ladder_crossing_count(ladder_nodes(patch)), 1L)
  expect_equal(ladder_crossing_count(
    ladder_nodes(ladder_patch_two_by_two(cross = FALSE))), 0L)

  # The other half of adversariality, and it was measured nowhere: the two species
  # must differ by whole factors in every per-species reduction parameter. A sum
  # that collapses them into one scalar, and a scatter that reaches the wrong
  # species, are both invisible while the two carry commensurate values.
  message(sprintf("  per-species reduction-parameter ratios: %s",
                  paste(sprintf("%.3f", ladder_species_ratio(patch)),
                        collapse = " ")))
  expect_true(ladder_species_separated(patch))

  # A stand's state is whatever the trajectory reached, so two of the table's
  # entries are measured there rather than required: the reserve band and the gate
  # slope. They are printed because that is the whole of their value -- a run stand
  # sits near half a relative reserve with a gate slope an order below the floor,
  # so every channel running through growth is damped relative to the constructed
  # patch the block Jacobian is formed on, and a check that never prints the number
  # asserts the regime instead of reporting it.
  stand <- ladder_stand_two_by_two()
  stand_report <- ladder_regime_report(stand, "stand")
  message("\nstand fixture regime:")
  for (i in seq_len(nrow(stand_report))) {
    message(sprintf("  %-58s %-5s %-7s %s", stand_report$assertion[[i]],
                    stand_report$ok[[i]],
                    if (stand_report$enforced[[i]]) "" else "measured",
                    stand_report$value[[i]]))
  }
  expect_true(all(stand_report$ok[stand_report$enforced]))

  # The floor a tolerance is measured against, rather than a literal. It is also
  # the cheapest form of the purity check: a forward quantity that does not come
  # back to its own value after the shared solver has been driven elsewhere is
  # carrying something between solves.
  floor <- ladder_forward_floor(patch)
  message(sprintf("  forward arithmetic floor: %.3e", floor))
  expect_lt(floor, 1e-12)
})

test_that("the census integrates over the coordinate the density is carried in", {
  # A gradient of a wrong number is worthless. A census is a quadrature of a
  # density, so its weights must be gaps in the coordinate the state is carried
  # on. Where that is the birth date, a grid built from heights integrates a
  # density in one variable against the spacing of another -- wrong in value on
  # any stand, and merely more visibly wrong on a crossed one.
  stand <- ladder_stand_two_by_two()
  ladder_require_regime(stand, "stand")

  model <- unname(stand_census(stand)["leaf_area"])
  over_birth_date <- ladder_census_leaf_area(stand, "birth_date")
  over_height <- ladder_census_leaf_area(stand, "height")

  # Non-vacuity: the two axes must disagree, or this check would pass on a stand
  # that cannot distinguish them.
  expect_gt(abs(over_birth_date - over_height) / abs(over_height), 1e-3)
  message(sprintf(
    "\n  census leaf_area: model %.8g | over birth date %.8g | over height %.8g",
    model, over_birth_date, over_height))

  expect_equal(model, over_birth_date, tolerance = 1e-8)
})

test_that("the census refuses a grid its own abscissa cannot order", {
  # Birth dates are strictly increasing by construction, so on this coordinate
  # the guard is an assertion and no sort is needed. What must not happen is a
  # census that silently integrates a tied or inverted grid.
  patch <- ladder_patch_two_by_two()
  nodes <- ladder_nodes(patch)
  for (i in unique(nodes$species)) {
    part <- nodes[nodes$species == i & !nodes$boundary, , drop = FALSE]
    expect_false(any(duplicated(part$birth_date)))
  }
  # And the boundary node is the youngest of all, since it is the one about to
  # be introduced. A reduction that placed it anywhere else would be integrating
  # over the wrong endpoint.
  for (i in unique(nodes$species)) {
    part <- nodes[nodes$species == i, , drop = FALSE]
    expect_equal(which.max(part$birth_date), which(part$boundary))
  }
})

test_that("the reference computes the model's own rates", {
  # The cheapest check in the corpus and the one everything above it rests on: a
  # reference whose VALUE disagrees with the model is a reference to a different
  # function, and every Jacobian taken from it is that function's Jacobian --
  # internally consistent, plausible, and refereeing nothing.
  #
  # It needs no derivative and no seed, so it is a floor item rather than part of
  # the harness, and it is not subsumed by any agreement between a tangent and a
  # sweep: both are built on the same rebound twin, so a twin that differs from
  # the model in value can have them agree with each other exactly.
  for (fixture in list(ladder_patch_one(),
                       ladder_patch_two_by_two(cross = FALSE),
                       ladder_patch_two_by_two(cross = TRUE))) {
    model <- fixture$ode_rates
    reference <- ladder_rhs_value_forward_tf24(fixture)
    scale <- pmax(abs(model), abs(reference))
    relative <- ifelse(scale > 0, abs(model - reference) / scale, 0)
    worst <- which.max(relative)
    names(relative) <- ladder_rate_names(fixture)
    ladder_report_margin(
      sprintf("reference value, worst at %s", names(relative)[[worst]]),
      max(relative), 1e3 * ladder_forward_floor(fixture))
  }
})

# The sweep computes the same thing however it is decomposed.

test_that("the regime holds at the states a trajectory passed through", {
  # A regime report taken at the terminal state says nothing about the states the
  # trajectory passed through, and it is exactly those states the trajectory rungs
  # measure on. Every margin taken along a trajectory carries this qualification, so
  # the qualification is measured rather than assumed.
  #
  # A violated assertion invalidates the run rather than failing it, which is the
  # same convention the per-fixture report uses: the fixture is then outside the
  # domain its checks are declared for, and that is a finding about the fixture.
  stand <- ladder_stand_trajectory(lifetime = 2)
  trajectory <- stand$store_trajectory()
  patch <- ladder_as_patch(stand)
  at <- unique(round(seq(2, length(trajectory), length.out = 8)))

  # Which assertions are properties of a STATE rather than of the fixture's
  # construction. Adversariality is built into a fixture once -- distinct heights,
  # distinct moistures, separated species -- and a trajectory is entitled to
  # violate it: a cohort one step after its introduction sits at exactly the
  # boundary node's height, so "no two heights equal" is false there by
  # construction and says nothing about the regime the checks need.
  #
  # What 08 asks to be evaluated at every solve is the physical regime, and that is
  # this list.
  per_state <- c("soil moisture strictly interior on every layer",
                 "net production strictly positive at every cohort",
                 "light at every read at least 100 times the floor")

  collected <- list()
  for (k in at) {
    patch$set_ode_state(trajectory[[k]]$state, trajectory[[k]]$time)
    # The second evaluation of the inflow condition, which a state load alone does
    # not reach and which every quantity a regime is read off depends on.
    invisible(patch$ode_rates)
    report <- ladder_regime_report(patch, "stand")
    report <- report[report$assertion %in% per_state, , drop = FALSE]
    report$state <- k
    bad <- report[!report$ok, , drop = FALSE]
    if (nrow(bad) > 0L) {
      message(sprintf("  state %3d leaves the regime: %s", k,
                      paste(sprintf("%s (%s)", bad$assertion, bad$value),
                            collapse = "; ")))
    }
    collected[[length(collected) + 1L]] <- report
  }

  # The whole point is the range over the trajectory rather than any one state, so
  # the assertions are collected over every sampled state and reported once.
  enforced <- do.call(rbind, collected)
  expect_equal(nrow(enforced), length(at) * length(per_state))
  failures <- enforced[!enforced$ok, , drop = FALSE]
  message(sprintf("\n  %d enforced assertions over %d sampled states, %d violated",
                  nrow(enforced), length(at), nrow(failures)))
  if (nrow(failures) > 0L) {
    skip(paste("the trajectory leaves the declared regime at some sampled state,",
               "so runs on it are invalid rather than failing:",
               paste(unique(failures$assertion), collapse = "; ")))
  }
  expect_equal(nrow(failures), 0L)
})

test_that("a cohort's rates are reproducible from its own boundary", {
  # The premise the whole cohort-granular decomposition rests on, and the one
  # nothing structural defends: a cohort's rates are a pure function of its own
  # state, the environment values it reads, and the traits. Every cohort of a
  # species points at one strategy and therefore one leaf, so anything that object
  # remembers between calls is a channel the reverse pass cannot see.
  #
  # The block re-runs one individual on a FRESH copy of the strategy, from that
  # node's own packed boundary. In isolation there is no predecessor, so a buffer
  # sized but not cleared, an exit that writes part of an operating point, and a
  # cache keyed on less than it depends on all show up here and nowhere else -- the
  # forward pass is order-deterministic, so each of them reproduces exactly and
  # every double-valued test passes.
  for (label in c("one cohort", "four nodes", "three of one species")) {
    patch <- switch(label,
                    "one cohort" = ladder_patch_one(),
                    "four nodes" = ladder_patch_two_by_two(),
                    "three of one species" = ladder_patch_permutable())
    table <- ladder_node_rate_table(patch)
    # The walk above reads `sp$nodes`; `sp$size` is the count the patch keeps.
    # Comparing them is what says the walk saw every node rather than most of them.
    expect_equal(length(table),
                 sum(vapply(patch$species, function(s) as.numeric(s$size), 0)))
    # A hundred times the fixture's own floor, not ten. The block evaluates on the
    # active scalar and the patch on the double one, so the two differ by however
    # the compiler contracts each; measured worst is ten machine units on the
    # four-node patch. A carried quantity is another plant's operating point and is
    # orders above either.
    bound <- 100 * ladder_forward_floor(patch)

    worst <- 0
    for (k in seq_along(table)) {
      own <- ladder_block_value_tf24(patch, k)[seq_along(table[[k]]$rates)]
      scale <- pmax(abs(table[[k]]$rates), .Machine$double.xmin)
      worst <- max(worst, max(abs(own - table[[k]]$rates) / scale))
    }
    ladder_report_margin(paste("rates from the node's own boundary,", label),
                         worst, bound)
  }
})

test_that("the same plants solved in the opposite order give the same rates", {
  # The permutation, which is the only check a reordering can fail and a re-run
  # cannot. The forward pass is order-deterministic, so a stale read reproduces
  # exactly; what it cannot survive is being asked to solve the same two plants in
  # the other order.
  #
  # The two grid points exchanged carry equal quadrature weights by construction,
  # so the reductions see the same sum in the reverse order. That leaves the field's
  # knot values agreeing in their last bits rather than exactly, which is the floor
  # this check is measured against rather than a tolerance anyone picked.
  a <- ladder_patch_permutable(FALSE)
  b <- ladder_patch_permutable(TRUE)
  ta <- ladder_node_rate_table(a)
  tb <- ladder_node_rate_table(b)

  floor <- ladder_permutation_floor(a, b)
  message(sprintf("\n  the permutation moves the field's knot values by %.2e", floor))

  moved <- 0L
  worst <- 0
  for (x in ta) {
    at <- which(vapply(tb, function(y) y$height == x$height, logical(1)))
    expect_length(at, 1L)
    y <- tb[[at]]
    if (y$position != x$position) moved <- moved + 1L
    scale <- pmax(abs(x$rates), .Machine$double.xmin)
    worst <- max(worst, max(abs(x$rates - y$rates) / scale))
  }
  # Non-vacuity: if no state changed position the two runs are one run and the
  # equality below holds for a reason that says nothing about a carried quantity.
  expect_gte(moved, 2L)
  message(sprintf("  %d of %d states solved at a different position", moved,
                  length(ta)))
  ladder_report_margin("rates under a permutation of the solve order", worst,
                       100 * floor)
})

test_that("a rejected step attempt is not a question about the gradient", {
  # This is where the ladder asked for the rejection count, and the requirement is
  # withdrawn rather than skipped. The reverse pass never sees a rejected attempt:
  # the adaptive pass resolves the schedule, one state and one step size are
  # recorded per ACCEPTED step, and the sweep runs over those. There is no
  # exclusion left to perform, so "the gradient is unchanged when rejections are
  # excluded" compares a thing with itself.
  #
  # The seam a rejected attempt does leave is in the FORWARD replay, not the
  # sweep: it moves patch state that is not ODE state -- the first-same-as-last
  # derivative carry, and anything cached on the patch -- so a pinned replay is
  # not bit-identical to the adaptive run it replays. That is what the trajectory
  # reference's own floor measures, and it is asserted where the reference is
  # used rather than here.
  shared <- ladder_shared("introductions")
  stand <- shared$stand
  result <- shared$gradient
  ladder_report_margin("the replay reaches this stand's own census",
                       ladder_replay_floor(stand, result$value), 1e-9)
})

test_that("a sweep that never ran is distinguishable from an insensitive stand", {
  # An empty segment list is not an insensitive stand. Both must not come back
  # as a vector of zeros.
  p <- ladder_parameters("fast")
  p$node_schedule_times <- list(numeric(0))
  empty <- ladder_run(p)
  expect_error(stand_gradient(empty))
})

test_that("the reverse pass refuses the coordinate it cannot transpose", {
  # On the height coordinate the abscissa is state, so the reduction transposes
  # would omit a weight derivative and the recorded step would omit the density
  # rate's compression term. The sweep is then the transpose of a function the
  # forward model is not evaluating, and nothing about the arithmetic complains.
  p <- ladder_parameters("fast")
  p$node_schedule_times <- list(c(0, 0.63))
  stand <- ladder_run(p, ctrl = Control(node_density_in_birth_date = FALSE))

  expect_error(stand_census_state_adjoint(stand), "birth-date")
  expect_error(stand_gradient(stand), "birth-date")
  # The census value is not a gradient and is answerable on either coordinate,
  # so refusing it too would refuse a defined answer.
  expect_silent(stand_census(stand))
})

# Every zero is attributable to a named cause.

test_that("every trait column resolves to one declared class", {
  # An exact zero is the signature of a missing accumulator and never of true
  # insensitivity, which makes it indistinguishable from an ecological finding.
  # So zeros are classified rather than tolerated, and the classification is two
  # failure modes and not one: an exact zero in the live class is a missing
  # accumulator, and a registered parameter reaching no equation comes back as
  # round-off, which reads as a gradient.
  shared <- ladder_shared("two_by_two")
  stand <- shared$stand
  result <- shared$gradient
  g <- result$gradient

  declared_zero <- ladder_zero_by_construction()
  band <- ladder_roundoff_band()
  bare <- ladder_bare_traits(colnames(g))

  classify <- function(column, name) {
    peak <- max(abs(column))
    declared <- name %in% declared_zero
    if (declared && peak == 0) return("zero by construction")
    # A declared zero that comes back non-zero is as much a finding as an
    # undeclared zero: the declaration names a reason, and a number means the
    # reason has stopped holding.
    if (declared) return("declared zero, now live")
    if (peak == 0) return("exact zero, unaccounted")
    if (peak < band[[2]]) return("round-off, reaches no equation")
    "live"
  }

  classes <- vapply(seq_along(bare),
                    function(j) classify(g[, j], bare[[j]]), character(1))
  names(classes) <- colnames(g)

  message("\ntrait column classes:")
  for (cls in unique(classes)) {
    message(sprintf("  %-30s %s", cls,
                    paste(names(classes)[classes == cls], collapse = " ")))
  }

  expect_length(names(classes)[classes == "declared zero, now live"], 0L)
  expect_length(names(classes)[classes == "exact zero, unaccounted"], 0L)
  expect_length(names(classes)[classes == "round-off, reaches no equation"], 0L)
  # Non-vacuity: a declaration that named nothing would make the class above
  # unfalsifiable, and every column being live would make this check pass with
  # the classification absent.
  expect_gt(sum(classes == "zero by construction"), 0L)
  expect_gt(sum(classes == "live"), 0L)
})

test_that("each declared zero is zero for the cause it is declared for", {
  # Naming a cause is only worth more than tolerating a zero if the cause is
  # measured. Each cause is a claim about which rows of the right-hand side the
  # parameter moves, and the trait Jacobian of one rate evaluation answers it
  # directly: a column of zeros for a parameter claimed to move one rate, or a
  # moved rate for one claimed to move none, is a wrong declaration.
  # One species, because the cause is a property of the parameter and not of
  # which species carries it, and the trait Jacobian costs one evaluation per
  # column: forty-four here against eighty-eight on the four-node fixture, for
  # the same claim. The per-species question belongs to rungs 3 and 4.
  patch <- ladder_patch_one()
  ladder_require_regime(patch, "patch")
  j <- ladder_rhs_trait_jacobian_forward_tf24(patch)
  columns <- ladder_trait_names_tf24(patch)
  bare <- ladder_bare_traits(columns)
  rates <- ladder_rate_names(patch)

  message("\ndeclared zeros, and the rates they move:")
  for (name in ladder_zero_by_construction()) {
    # Both species' columns, so a cause that holds for one species and not the
    # other is a failure rather than a pass on the first match.
    at <- which(bare == name)
    expect_length(at, length(patch$species))
    moved <- rates[apply(abs(j[, at, drop = FALSE]) > 0, 1, any)]
    message(sprintf("  %-16s %-46s %s", name, ladder_zero_cause(name),
                    paste(unique(moved), collapse = " ")))
    if (name %in% ladder_zero_at_an_interior_optimum()) {
      # The dry bound of a feasible interval the operating point is inside. It
      # moves nothing at all, and that is complementary slackness rather than a
      # missing row -- the same parameter carries the whole row at a pin.
      expect_length(moved, 0L)
    } else {
      # The two reproductive accumulators and nothing else. A third rate would
      # make the census's silence about the column wrong, and neither of these
      # two is read by any metric or by any equation -- which is why the same
      # column would be live on a fitness functional.
      expect_setequal(unique(moved), ladder_reproductive_rates())
    }
  }
})

test_that("an unknown parameter refuses by name rather than returning a number", {
  # A zero and an absence are different at this boundary. An unknown parameter
  # must be refused; a registered one that reaches nothing is the other class
  # above and is the one that reads as an answer.
  stand <- ladder_stand_two_by_two()
  expect_error(stand_gradient(stand, traits = "not_a_parameter"), "Unknown trait")
  expect_error(stand_gradient(stand, metrics = "not_a_metric"),
               "Unknown census metric")
})

test_that("gradient columns are named per species", {
  # Concatenating each species' parameter names with no species prefix yields one
  # name repeated per species. Character indexing then resolves each to its first
  # match, so a multi-species gradient silently returns species one's column for
  # every named parameter, and unknown-parameter validation cannot see it. A
  # single-species suite never detects this.
  stand <- ladder_stand_two_by_two()
  columns <- census_trait_names_tf24(stand)
  expect_equal(length(columns),
               2L * length(unique(sub("^[0-9]+[.]", "", columns))))
  expect_false(any(duplicated(columns)))
})
