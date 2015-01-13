if (interactive()) {
  devtools::load_all("../../")
  library(testthat)
  source("helper-tree2.R")
}

context("Plant")

test_that("", {
  cmp <- make_reference_plant()
  s <- Strategy()
  p <- Plant(s)

  expect_that(p$strategy, is_identical_to(s))

  ## Set the height to something (here 10)
  h0 <- 10
  p$height <- h0

  p_size <- p$vars_size

  expect_that(p_size[["height"]],
              is_identical_to(h0))
  expect_that(p_size[["leaf_area"]],
              equals(cmp$LeafArea(h0)))
  expect_that(p_size[["mass_leaf"]],
              equals(cmp$LeafMass(cmp$traits$lma, cmp$LeafArea(h0))))
  expect_that(p_size[["mass_sapwood"]],
              equals(cmp$SapwoodMass(cmp$traits$rho, cmp$LeafArea(h0), h0)))
  expect_that(p_size[["mass_bark"]],
              equals(cmp$BarkMass(cmp$traits$rho, cmp$LeafArea(h0), h0)))
  expect_that(p_size[["mass_root"]],
              equals(cmp$RootMass(cmp$LeafArea(h0))))
  expect_that(p_size[["mass_live"]],
              equals(cmp$LiveMass(cmp$traits, cmp$LeafArea(h0))))
  expect_that(p_size[["area_bark"]],
              equals(cmp$bark_area(h0)))
  expect_that(p_size[["area_sapwood"]],
              equals(cmp$sapwood_area(h0)))
  expect_that(p_size[["area_basal"]],
              equals(cmp$basal_area(h0)))

  expect_that(p$height,    is_identical_to(p_size[["height"]]))
  expect_that(p$leaf_area, is_identical_to(p_size[["leaf_area"]]))

  ## Heartwood function
  ## TODO: Check with Daniel here -- might need updating.
  ## Expect zero unless it has been set otherwise
  expect_that(p_size[["area_heartwood"]],
              is_identical_to(0.0))
  HA0 <- 1e-3
  p$heartwood_area <- HA0
  p_size <- p$vars_size
  expect_that(p_size[["area_heartwood"]],
              is_identical_to(HA0))
  expect_that(p_size[["area_basal"]],
              equals(cmp$basal_area(h0) + HA0))
  # set heartwood back at zero for subsequent tests
  p$heartwood_area <- 0
})
