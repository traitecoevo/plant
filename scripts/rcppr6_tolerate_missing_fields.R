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
## A Control saved before node_density_coordinate replaced the boolean
## node_density_in_birth_date (#516) carries the boolean. Tolerating the missing
## field alone would read it as "auto", which for a size-only model is the
## height coordinate whatever the boolean chose, so read the boolean as the
## explicit coordinate it named.
legacy <- c(
  "  // A Control saved before node_density_coordinate replaced the boolean",
  "  // node_density_in_birth_date carries the boolean: read it as the explicit",
  "  // coordinate it named, not as \"auto\".",
  "  if (!xl.containsElementNamed(\"node_density_coordinate\") &&",
  "      xl.containsElementNamed(\"node_density_in_birth_date\")) {",
  "    ret.node_density_coordinate =",
  "      Rcpp::as<bool>(xl[\"node_density_in_birth_date\"]) ? \"birth_date\" : \"height\";",
  "  }")
control <- grep("^template <> inline plant::Control as\\(SEXP x\\) \\{$", lines)
if (length(control) != 1L) {
  stop("expected one Control converter in ", path)
}
if (!any(grepl("node_density_in_birth_date\"]) ?", lines, fixed = TRUE))) {
  i <- control
  while (!grepl("^\\s*return ret;$", lines[i])) {
    i <- i + 1L
  }
  lines <- append(lines, legacy, after = i - 1L)
}
writeLines(lines, path)
message("tolerated missing fields in ", length(starts), " converters (",
        n, " fields)")
