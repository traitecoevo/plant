# Extracted from test-strategy-tf24.R:97

# test -------------------------------------------------------------------------
s <- TF24_Strategy()
p <- TF24_Individual(s)
expect_equal(p$aux_size, 13)
