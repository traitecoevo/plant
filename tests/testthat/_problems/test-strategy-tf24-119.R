# Extracted from test-strategy-tf24.R:119

# test -------------------------------------------------------------------------
s <- TF24_Strategy()
p <- TF24_Individual(s)
expect_equal(p$aux_size, 13)
expect_equal(length(p$internals$auxs), 13)
expect_equal(p$aux_names, c(
    "competition_effect",
    "height_inverse",
    "net_mass_production_dt",
    "root_mass",
    "opt_psi_stem",
    "opt_root_psi",
    "transpiration",
    "E_up_",
    "profit",
    "shadow_cost",
    "stom_cond_CO2",
    "assimilation",
    "leaf_marginal_return"
  ))
s <- TF24_Strategy(collect_all_auxiliary=TRUE)
expect_true(s$collect_all_auxiliary)
p <- TF24_Individual(s)
expect_equal(p$aux_size, 14)
expect_equal(length(p$internals$auxs), 14)
