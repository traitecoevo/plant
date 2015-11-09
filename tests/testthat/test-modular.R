context("Modular")

test_that("Construction", {

  ## This is a *minimal* set of tests that checks that it is possible
  ## to create the full set of different object types (Plant, Species,
  ## Patch etc) with different underlying Strategy types.  It doesn't
  ## actually try to run them though, so do that elsewhere.

  strategy_types <- get_list_of_strategy_types()

  for (x in names(strategy_types)) {
    s <- strategy_types[[x]]()
    expect_that(s, is_a(paste0(x, "_Strategy")))

    p <- Plant(x)(s)

    expect_that(p, is_a("Plant"))
    expect_that(p, is_a(sprintf("Plant<%s>", x)))
    expect_that(class(p$strategy), equals(class(s)))
    expect_that(p$strategy, equals(s))

    pp <- PlantPlus(x)(s)
    expect_that(pp, is_a("PlantPlus"))
    expect_that(pp, is_a(sprintf("PlantPlus<%s>", x)))
    expect_that(class(pp$strategy), equals(class(s)))
    expect_that(pp$strategy, equals(s))

    coh <- Cohort(x)(s)

    expect_that(coh, is_a("Cohort"))
    expect_that(coh, is_a(sprintf("Cohort<%s>", x)))
    expect_that(class(coh$plant), equals(class(p)))

    sp <- Species(x)(s)

    expect_that(sp, is_a("Species"))
    expect_that(sp, is_a(sprintf("Species<%s>", x)))
    expect_that(class(sp$seed), equals(class(coh)))

    par <- Parameters(x)(strategies=list(s))
    expect_that(par, is_a("Parameters"))
    expect_that(par, is_a(sprintf("Parameters<%s>", x)))
    expect_that(par$strategies[[1]], equals(s))

    pat <- Patch(x)(par)
    expect_that(pat, is_a("Patch"))
    expect_that(pat, is_a(sprintf("Patch<%s>", x)))
    expect_that(class(pat$species[[1]]), equals(class(sp)))

    ebt <- SCM(x)(par)
    expect_that(ebt, is_a("SCM"))
    expect_that(ebt, is_a(sprintf("SCM<%s>", x)))
    expect_that(class(ebt$patch), equals(class(pat)))

    ## Stochastic model:
    s_sp <- StochasticSpecies(x)(s)
    expect_that(s_sp, is_a("StochasticSpecies"))
    expect_that(s_sp, is_a(sprintf("StochasticSpecies<%s>", x)))
    expect_that(class(s_sp$seed), equals(class(p)))

    s_pat <- StochasticPatch(x)(par)
    expect_that(s_pat, is_a("StochasticPatch"))
    expect_that(s_pat, is_a(sprintf("StochasticPatch<%s>", x)))
    expect_that(class(s_pat$species[[1]]), equals(class(s_sp)))

    s_pr <- StochasticPatchRunner(x)(par)
    expect_that(s_pr, is_a("StochasticPatchRunner"))
    expect_that(s_pr, is_a(sprintf("StochasticPatchRunner<%s>", x)))
    expect_that(class(s_pr$patch), equals(class(s_pat)))
  }
})
