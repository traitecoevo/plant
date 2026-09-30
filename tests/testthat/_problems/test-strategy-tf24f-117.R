# Extracted from test-strategy-tf24f.R:117

# test -------------------------------------------------------------------------
s <- TF24f_Strategy()
p <- TF24f_Individual(s)
expect_equal(p$aux_size, 13)
expect_equal(length(p$internals$auxs), 13)
