
source("tests/testthat/helper-plant.R")
x <- "FF16r"
s <- strategy_types[[x]]()
p <- Parameters(x, paste0(x, "_Env"))(hyperpar=FF16_hyperpar)

expect_is(p, sprintf("Parameters<%s,%s_Env>", x, x))
p$cohort_schedule_max_time
cohort_schedule_max_time_default(p)
cohort_schedule_max_time_default__Parameters___FF16__FF16_Env(p)
