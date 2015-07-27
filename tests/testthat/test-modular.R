context("Modular")

test_that("Construction", {
  ## Do this in a generic way so that we can add as many additional
  ## classes in the future and test in one place.
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
  }
})
