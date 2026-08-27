# Trait names are checked, not silently ignored (#636).
#
# `generate_strategy()` assigns trait columns into the strategy's nested `pars`.
# `pars[["typo"]] <- x` on a list APPENDS, so before #636 an unknown name was not
# dropped: it became a junk element that the generated `Rcpp::as<*_Pars>` never
# looks at (`RcppR6_post.hpp`: "No current support for a hook"). The strategy then
# built, solved and reported at the DEFAULT value of the parameter the caller
# believed they were varying, with nothing in the output to say so -- a trait sweep
# running its whole length at one point, or a calibration fitting a model it never
# perturbed.
#
# The valid set is deliberately wider than `names(pars)`: a hyperpar may CONSUME a
# trait it never assigns, and it declares those itself via the `input_traits`
# attribute. FF16's `narea` is the case -- it derives `a_p1`, `a_p2` and `r_l` and
# is not a parameter of any model.

test_that("an unknown trait name is refused, not ignored", {
  p0 <- scm_base_parameters("FF16")
  # The regression: 'lmaa' used to be accepted and lma stayed at its default.
  expect_error(generate_strategy(p0, trait_matrix(0.5, "lmaa")),
               "Unknown trait name for FF16_Strategy")
  expect_error(generate_strategy(p0, trait_matrix(0.5, "lmaa")),
               "did you mean lma")
})

test_that("the error names the model and every offender", {
  p0 <- scm_base_parameters("FF16")
  m <- trait_matrix(matrix(c(0.5, 600, 1), nrow = 1),
                    c("lma", "rhoo", "nonsense_xyz"))
  err <- tryCatch(generate_strategy(p0, m), error = conditionMessage)
  # Plural, both offenders, and NOT the one valid column.
  expect_match(err, "Unknown trait names")
  expect_match(err, "rhoo")
  expect_match(err, "nonsense_xyz")
  expect_false(grepl("\\blma\\b \\(did", err))
  expect_match(err, "FF16_Strategy")
  # A name with no near match gets no misleading suggestion.
  expect_false(grepl("nonsense_xyz \\(did you mean", err))
})

test_that("the hyperpar's own arguments are refused, and said to be arguments", {
  # ⚠️ THE MISTAKE WORTH CATCHING BY NAME. `B_kl1` and friends parameterise the
  # trade-offs once, when the hyperpar is BUILT; they are not traits.
  # test-ode-individual-runner.R passed nine of them for years to no effect.
  p0 <- scm_base_parameters("FF16")
  expect_true("B_kl1" %in% names(formals(make_FF16_hyperpar)))
  expect_false("B_kl1" %in% names(p0$strategy_default$pars))
  err <- tryCatch(generate_strategy(p0, trait_matrix(0.5, "B_kl1")),
                  error = conditionMessage)
  expect_match(err, "make_FF16_hyperpar\\(\\) instead", fixed = FALSE)
})

test_that("every model refuses an unknown name", {
  for (type in c("FF16", "K93", "TF24", "TF24f")) {
    p0 <- scm_base_parameters(type)
    expect_error(generate_strategy(p0, trait_matrix(1, "definitely_not_a_trait")),
                 sprintf("Unknown trait name for %s_Strategy", type),
                 info = type)
  }
})

test_that("a legitimate trait sweep is unaffected", {
  p0 <- scm_base_parameters("FF16")
  lma <- c(0.05, 0.1, 0.5)
  ss <- generate_strategy(p0, trait_matrix(lma, "lma"))
  expect_equal(vapply(ss, function(s) s$pars$lma, numeric(1)), lma)
  # And nothing spurious is left on `pars`.
  expect_setequal(names(ss[[1]]$pars), names(p0$strategy_default$pars))
})

test_that("a trait the hyperpar consumes but never assigns is accepted", {
  # ⚠️ WITHOUT THE `input_traits` DECLARATION THIS IS THE FALSE POSITIVE the guard
  # would introduce: `narea` is not a parameter of FF16, so a `pars`-only check
  # rejects a sweep that works and matters.
  p0 <- scm_base_parameters("FF16")
  expect_false("narea" %in% names(p0$strategy_default$pars))
  expect_identical(attr(make_FF16_hyperpar(), "input_traits"), "narea")

  narea <- c(0.001, 0.01)
  ss <- generate_strategy(p0, trait_matrix(narea, "narea"))
  # It did its work: narea drives a_p1, a_p2 and r_l.
  expect_true(all(diff(vapply(ss, function(s) s$pars$a_p1, numeric(1))) > 0))
  expect_true(all(diff(vapply(ss, function(s) s$pars$r_l,  numeric(1))) > 0))
  # ...without being appended to pars as a junk element.
  expect_setequal(names(ss[[1]]$pars), names(p0$strategy_default$pars))
})

test_that("pars and consumed traits can be mixed in one matrix", {
  # The positional hazard: `xi` is unnamed and aligned with the column order, so
  # assigning by name would silently write NAs. A mixed matrix, where the assigned
  # columns are not a prefix of the trait names, is what catches that.
  p0 <- scm_base_parameters("FF16")
  m <- trait_matrix(matrix(c(0.005, 0.5, 700), nrow = 1),
                    c("narea", "lma", "rho"))
  s <- generate_strategy(p0, m)[[1]]
  expect_equal(s$pars$lma, 0.5)
  expect_equal(s$pars$rho, 700)
  expect_true(is.finite(s$pars$a_p1))
  expect_false(anyNA(unlist(unclass(s$pars))))
})

test_that("a caller-supplied hyperpar with no declaration gets pars alone", {
  # We cannot know what an arbitrary closure reads, so the honest default is the
  # narrow set. Documented rather than merely true.
  p0 <- scm_base_parameters("FF16")
  bare <- function(m, s, filter = TRUE) m
  expect_null(attr(bare, "input_traits", exact = TRUE))
  expect_silent(generate_strategy(p0, trait_matrix(0.5, "lma"), hyperpar = bare))
  expect_error(generate_strategy(p0, trait_matrix(0.005, "narea"),
                                 hyperpar = bare),
               "Unknown trait name")
})

test_that("add_strategies and add_mutant refuse it too", {
  # They route through generate_strategy, so this is about the routes a user
  # actually takes rather than about a second check.
  p0 <- scm_base_parameters("FF16")
  expect_error(add_strategies(p0, trait_matrix(0.5, "lmaa"), birth_rate = 1),
               "Unknown trait name")
  expect_error(add_mutant(p0, trait_matrix(0.5, "lmaa"), birth_rate = 1),
               "Unknown trait name")
})
