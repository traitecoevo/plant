context("Modular")

test_that("Construction", {

  ## This is a *minimal* set of tests that checks that it is possible
  ## to create the full set of different object types (Plant, Species,
  ## Patch etc) with different underlying Strategy types.  It doesn't
  ## actually try to run them though, so do that elsewhere.

  cl <- list(FFW16=FFW16_Strategy,
             FFdev=FFdev_Strategy)

  for (x in names(cl)) {
    s <- cl[[x]]()
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

    ebt <- EBT(x)(par)
    expect_that(ebt, is_a("EBT"))
    expect_that(ebt, is_a(sprintf("EBT<%s>", x)))
    expect_that(class(ebt$patch), equals(class(pat)))
  }
})
