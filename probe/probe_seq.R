## Rscript probe_seq.R <body.R> <label>: the test body as written, then the
## invader gap; and again with each earlier call on `substepped` dropped.
suppressPackageStartupMessages({ library(plant); library(testthat) })
body <- readLines(commandArgs(TRUE)[1])
label <- commandArgs(TRUE)[2]
cut <- grep("^  ratios <- ", body)
pre <- body[seq_len(cut - 1)]
gap <- function(lines, what) {
  e <- new.env(parent = asNamespace("plant"))
  eval(parse(text = lines), e)
  r0 <- e$substepped$net_reproduction_ratios
  e$substepped$run_mutant(e$p)
  r1 <- e$substepped$net_reproduction_ratios
  cat(sprintf("%-8s %-26s stand %.17g invader %.17g rel %.3e\n", label, what,
              r0, r1, abs(r1 - r0) / abs(r0)))
}
drop <- function(pattern) {
  hit <- grep(pattern, pre)
  out <- pre
  out[hit] <- paste0("# ", out[hit])
  out
}
gap(pre, "as written")
gap(drop("^  g <- stand_gradient|expect_identical\\(names|expect_equal\\(unname\\(g\\$|^ +tolerance = 1e-2"),
    "no stand_gradient")
gap(drop("^  tangent <- ladder|^ +as.numeric\\(traits|expect_equal\\(tangent|tangent\\$tangent"),
    "no tangent")
gap(drop("^  differenced <-|expect_equal\\(unname\\(g\\$gradient\\[, \"1.lma\"\\]\\), unname\\(differenced|^ +tolerance = 1e-2"),
    "no differences")
