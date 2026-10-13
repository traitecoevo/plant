## Post-process RcppR6's generated converters so that a Control or Parameters
## object saved before a field was added still loads, with that field at its
## default. Run by `make RcppR6` after RcppR6::RcppR6(); safe to run twice.
##
## RcppR6 converts a list class by default-constructing the C++ object and then
## assigning every field from the R list, so a list missing a field fails with
## "Index out of bounds". New fields are opt-in and inert at their default, so
## an object that predates one means the default; reading it as such keeps old
## saved objects loadable (#661 review). Limited to the two classes saved and
## reloaded in practice, so a hand-built list of any other class still fails
## loudly on a missing field.

path <- "inst/include/plant/RcppR6_post.hpp"
lines <- readLines(path)

starts <- grep("^template <> inline plant::(Control|Parameters<.*>) as\\(SEXP x\\) \\{$",
               lines)
if (length(starts) == 0L) {
  stop("no Control or Parameters converter found in ", path)
}
field <- "^(\\s*)(ret\\.\\w+ = Rcpp::as<.*>\\(xl\\[\"(\\w+)\"\\]\\);)$"
n <- 0L
for (s in starts) {
  i <- s
  while (!grepl("^\\s*return ret;$", lines[i])) {
    if (grepl(field, lines[i])) {
      lines[i] <- sub(field, "\\1if (xl.containsElementNamed(\"\\3\")) \\2",
                      lines[i])
      n <- n + 1L
    }
    i <- i + 1L
  }
}
writeLines(lines, path)
message("tolerated missing fields in ", length(starts), " converters (",
        n, " fields)")
