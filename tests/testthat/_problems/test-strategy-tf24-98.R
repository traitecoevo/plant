# Extracted from test-strategy-tf24.R:98

# test -------------------------------------------------------------------------
s <- TF24_Strategy()
p <- TF24_Individual(s)
expect_equal(p$aux_size, 13)
expect_equal(length(p$internals$auxs), 13)
