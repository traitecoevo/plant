# Extracted from test-strategy-tf24f.R:116

# test -------------------------------------------------------------------------
s <- TF24f_Strategy()
p <- TF24f_Individual(s)
expect_equal(p$aux_size, 13)
