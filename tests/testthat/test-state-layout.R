# HEIGHT_INDEX, MORTALITY_INDEX and FECUNDITY_INDEX are a claim about the first
# three state slots of every model, made by around fifty readers across the
# package. Every model happens to satisfy it, so a test that only builds models
# cannot tell a live check from a dead one -- it has to be handed a layout that
# breaks the claim.

layout_probe <- function() {
  testthat::skip_if_not_installed("Rcpp")
  odelia_inc <- file.path(find.package("odelia"), "include")
  odelia_so <- Sys.glob(file.path(find.package("odelia"), "libs", "*odelia.so"))
  # Installed, then the source tree, because `load_all` leaves the library where
  # it was built rather than under libs/.
  plant_so <- c(Sys.glob(file.path(find.package("plant"), "libs", "*plant.so")),
                testthat::test_path("..", "..", "src", "plant.so"))
  plant_so <- plant_so[file.exists(plant_so)]
  testthat::skip_if(!dir.exists(odelia_inc), "odelia headers not found.")
  testthat::skip_if(!length(odelia_so) || !file.exists(odelia_so[[1]]),
                    "odelia shared library not found.")
  testthat::skip_if(!length(plant_so),
                    "plant shared library not found; util::stop lives in it.")
  inc <- c(system.file("include", package = "plant"), odelia_inc,
           file.path(find.package("BH"), "include"))
  inc <- inc[nzchar(inc)]

  # `// [[Rcpp::plugins(cpp20)]]` cannot go in PKG_CPPFLAGS: R places those
  # before its own -std= and wins, and every concept in the headers then reads
  # as a syntax error.
  withr::local_envvar(
    PKG_CPPFLAGS = paste(paste0("-I", shQuote(inc)), collapse = " "),
    PKG_LIBS = shQuote(normalizePath(plant_so[[1]])))
  res <- tryCatch({
    Rcpp::sourceCpp(code = '
      // [[Rcpp::plugins(cpp20)]]
      #include <Rcpp.h>
      #include <map>
      #include <string>
      #include <plant/strategy.h>

      // [[Rcpp::export]]
      std::string layout_says(int height_at, bool name_it) {
        std::map<std::string, int> m{{"mortality", 1}, {"fecundity", 2}};
        if (name_it) m["height"] = height_at;
        try {
          plant::check_state_layout(m, m.size(), "PROBE");
        } catch (const std::exception& e) {
          return std::string(e.what());
        }
        return "accepted";
      }

      // [[Rcpp::export]]
      std::string size_says(int claimed) {
        std::map<std::string, int> m{
            {"height", 0}, {"mortality", 1}, {"fecundity", 2}};
        try {
          plant::check_state_layout(m, static_cast<size_t>(claimed), "PROBE");
        } catch (const std::exception& e) {
          return std::string(e.what());
        }
        return "accepted";
      }', verbose = FALSE)
    NULL
  }, error = function(e) e)
  if (inherits(res, "error")) {
    testthat::skip(paste("the layout probe did not build:",
                         conditionMessage(res)))
  }
}

test_that("a model whose state order breaks the index claim is refused", {
  layout_probe()

  # The order every model does use.
  expect_equal(layout_says(0L, TRUE), "accepted")

  # Height named second. Unchecked, every generic reader of a height would read
  # a mortality instead -- a small positive number, so finite and plausible
  # wherever it lands.
  expect_match(layout_says(1L, TRUE), "names `height` state 1, against the 0")

  # And a model that names no height at all, which reads as position zero
  # holding whatever that model put first.
  expect_match(layout_says(0L, FALSE), "names no `height` state")
})

# state_size() is a literal on every model while the index map is built from
# state_names(), so the two can disagree. Unchecked, Internals is sized by the
# literal and a name past that length resolves to an index no vector holds.
test_that("a model whose state_size disagrees with its names is refused", {
  layout_probe()

  expect_equal(size_says(3L), "accepted")
  expect_match(size_says(2L), "names 3 states against the 2")
  expect_match(size_says(4L), "names 3 states against the 4")
})
