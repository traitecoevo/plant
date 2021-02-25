context("Modular")

test_that("Construction", {

  ## This is a *minimal* set of tests that checks that it is possible
  ## to create the full set of different object types (Individual, Species,
  ## Patch etc) with different underlying Strategy types.  It doesn't
  ## actually try to run them though, so do that elsewhere.

  strategy_types <- get_list_of_strategy_types()
  environment_types <- get_list_of_environment_types()

  for (x in names(strategy_types)) {
    s <- strategy_types[[x]]()
    e <- environment_types[[x]]
    expect_is(s, paste0(x, "_Strategy"))

    p <- Individual(x, e)(s)

    expect_is(p, "Individual")
    expect_is(p, sprintf("Individual<%s,%s>", x, e))
    expect_equal(class(p$strategy), class(s))
    expect_equal(p$strategy, s)

    coh <- Cohort(x, e)(s)

    expect_is(coh, "Cohort")
    expect_is(coh, sprintf("Cohort<%s,%s>", x, e))
    expect_equal(class(coh$plant), class(p))

    sp <- Species(x, e)(s)

    expect_is(sp, "Species")
    expect_is(sp, sprintf("Species<%s,%s>", x, e))
    expect_equal(class(sp$seed), class(coh))

    par <- Parameters(x, e)(strategies=list(s))
    expect_is(par, "Parameters")
    expect_is(par, sprintf("Parameters<%s,%s>", x, e))
    expect_equal(par$strategies[[1]], s)

    pat <- Patch(x, e)(par)
    expect_is(pat, "Patch")
    expect_is(pat, sprintf("Patch<%s,%s>", x, e))
    expect_equal(class(pat$species[[1]]), class(sp))

    scm <- SCM(x, e)(par)
    expect_is(scm, "SCM")
    expect_is(scm, sprintf("SCM<%s,%s>", x, e))
    expect_equal(class(scm$patch), class(pat))

    ## Stochastic model:
    s_sp <- StochasticSpecies(x, e)(s)
    expect_is(s_sp, "StochasticSpecies")
    expect_is(s_sp, sprintf("StochasticSpecies<%s,%s>", x, e))
    expect_equal(class(s_sp$seed), class(p))

    s_pat <- StochasticPatch(x, e)(par)
    expect_is(s_pat, "StochasticPatch")
    expect_is(s_pat, sprintf("StochasticPatch<%s,%s>", x, e))
    expect_equal(class(s_pat$species[[1]]), class(s_sp))

    s_pr <- StochasticPatchRunner(x, e)(par)
    expect_is(s_pr, "StochasticPatchRunner")
    expect_is(s_pr, sprintf("StochasticPatchRunner<%s,%s>", x, e))
    expect_equal(class(s_pr$patch), class(s_pat))
  }
})
