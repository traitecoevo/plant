context("SCM-general")


test_that("Basic functions", {
  
  l <- Leaf()
  
  K_s = 1
  h_v = 2
  h = 3
  
  expect_equal(l$calc_k_l_max(K_s, h_v, h), 2 / 3)
  
}
)
