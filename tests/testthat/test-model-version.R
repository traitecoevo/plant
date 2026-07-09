# Scientific model versioning (see R/strategy_support.R, src/strategy_version.cpp,
# and the `scientific_version` constant in inst/include/plant/models/*_strategy.h).

models <- c("FF16", "K93", "TF24", "TF24f")

test_that("model_version() returns a positive integer for every model", {
  for (type in models) {
    v <- model_version(type)
    expect_type(v, "integer")
    expect_length(v, 1L)
    expect_true(v >= 1L)
  }
})

test_that("model_id() is 'Model@vN' and agrees with model_version()", {
  for (type in models) {
    expect_identical(model_id(type),
                     sprintf("%s@v%d", type, model_version(type)))
  }
})

test_that("model_version() rejects an unknown model", {
  expect_error(model_version("NOPE"), "Unknown type")
})

# Drift guard -----------------------------------------------------------------
#
# Snapshot the default *scientific surface* of each model: its model_id plus the
# default strategy (biological `pars`, numerical `control`, and birth-rate
# defaults). If a default parameter changes but `scientific_version` is NOT
# bumped, the model_id stays the same while the snapshot content differs, so
# this test fails -- prompting the developer either to bump the version (in the
# model header) or to accept the change with testthat::snapshot_accept().
# Pure equation changes that move no default are not caught here and rely on
# reviewer discipline; see the bump convention in the header comments.
#
# Each default is pre-formatted to a full-precision string so the snapshot is
# floating-point exact, portable across R versions, and human-readable (unlike
# style = "serialize", which is R-version fragile and skipped on CRAN).

# Flatten a list of scalar defaults to a named character vector, formatting
# doubles with 17 significant digits so no precision is lost.
flatten_defaults <- function(lst, prefix) {
  out <- character(0)
  for (nm in names(lst)) {
    v <- lst[[nm]]
    key <- if (nzchar(prefix)) paste0(prefix, ".", nm) else nm
    if (length(v) == 0L) {
      out[key] <- ""
      next
    }
    vals <- vapply(v, function(x) {
      if (is.character(x)) x
      else if (is.logical(x)) as.character(x)
      else sprintf("%.17g", x)
    }, character(1L))
    names(vals) <- if (length(vals) > 1L) {
      paste0(key, "[", seq_along(vals), "]")
    } else {
      key
    }
    out <- c(out, vals)
  }
  out
}

scientific_surface <- function(type) {
  s <- as.list(get(paste0(type, "_Strategy"))())
  extras <- list(
    birth_rate_x = s$birth_rate_x,
    birth_rate_y = s$birth_rate_y,
    is_variable_birth_rate = s$is_variable_birth_rate
  )
  out <- c(
    model_id = model_id(type),
    flatten_defaults(s$pars, "pars"),
    flatten_defaults(s$control, "control"),
    flatten_defaults(extras, "")
  )
  out[order(names(out))]
}

for (type in models) {
  test_that(
    sprintf("scientific surface of %s matches its declared version", type), {
      expect_snapshot_value(scientific_surface(type), style = "json2")
    })
}
