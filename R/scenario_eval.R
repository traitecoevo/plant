##' Scenario evaluation framework for the TF24 / TF24f hydraulic models.
##'
##' These functions turn a table of qualitatively-described model scenarios
##' (e.g. \code{inst/scenarios/model_scenarios_hydraulic.csv}, added in PR #555)
##' into concrete TF24 parameterisations, run each through the SCM, classify the
##' run as a success or a failure, and score the observed outcome against the
##' expected outcome recorded in the table. The result is a provenance-stamped
##' scorecard that can be re-run across branches as a gateway check: many
##' scenarios are expected to fail on the current model and become targets to
##' fix as new features (e.g. NSC storage, #554) land.
##'
##' The qualitative-to-quantitative translation is data-driven: every
##' \code{High}/\code{Low}/descriptor cell is looked up in an editable mapping
##' table (\code{inst/scenarios/scenario_mapping.csv}) so the biology can be
##' recalibrated without touching code.
##'
##' @name scenario_eval
##' @family scenario_eval
NULL

## Path to a packaged scenario data file. Falls back to the source tree when
## the package is loaded via devtools/pkgload (system.file() then returns "").
scenario_file <- function(name) {
  path <- system.file("scenarios", name, package = "plant")
  if (!nzchar(path)) {
    path <- file.path("inst", "scenarios", name)
  }
  path
}

##' @param path Path to the CSV file.
##' @return \code{read_scenario_table} returns a tibble of scenarios with an
##'   added \code{scenario_id} and \code{is_duplicate} flag; the raw descriptor
##'   columns are kept verbatim.
##' @rdname scenario_eval
##' @export
read_scenario_table <- function(
    path = scenario_file("model_scenarios_hydraulic.csv")) {
  ## fileEncoding = "UTF-8-BOM" transparently strips a leading BOM if present.
  raw <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE,
                         fileEncoding = "UTF-8-BOM")
  descriptor_cols <- setdiff(names(raw), "Scenario")
  tbl <- tibble::as_tibble(raw)
  tbl$scenario_id <- sprintf("S%02d", seq_len(nrow(tbl)))
  tbl$is_duplicate <- duplicated(raw[descriptor_cols])
  ## Put the id first for readability.
  tbl[c("scenario_id", setdiff(names(tbl), "scenario_id"))]
}

##' @rdname scenario_eval
##' @export
read_scenario_mapping <- function(
    path = scenario_file("scenario_mapping.csv")) {
  raw <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE,
                         fileEncoding = "UTF-8-BOM")
  tibble::as_tibble(raw)
}

##' @param row A one-row data frame / tibble from \code{read_scenario_table}.
##' @param mapping A mapping tibble from \code{read_scenario_mapping}.
##' @return \code{scenario_to_config} returns a list with \code{traits} (named
##'   numeric), \code{env} (named list of Environment fields), \code{driver}
##'   (named list of extrinsic-driver settings) and \code{expected}
##'   (\code{"failure"} or \code{"success"}).
##' @rdname scenario_eval
##' @export
scenario_to_config <- function(row, mapping) {
  row <- as.list(row)
  traits <- list()
  env <- list()
  driver <- list()

  for (col in unique(mapping$csv_column)) {
    if (!col %in% names(row)) {
      next
    }
    level <- as.character(row[[col]])
    hit <- mapping[mapping$csv_column == col & mapping$level == level, ]
    if (nrow(hit) == 0L) {
      stop(sprintf(
        "No mapping for column '%s' level '%s' (scenario '%s'). Add a row to scenario_mapping.csv.",
        col, level, row$Scenario), call. = FALSE)
    }
    if (nrow(hit) > 1L) {
      stop(sprintf("Ambiguous mapping for column '%s' level '%s' (%d matches).",
                   col, level, nrow(hit)), call. = FALSE)
    }
    value <- as.numeric(hit$value)
    switch(hit$target,
           trait  = traits[[hit$param]]  <- value,
           env    = env[[hit$param]]     <- value,
           driver = driver[[hit$param]]  <- value,
           stop(sprintf("Unknown target '%s' in scenario_mapping.csv.", hit$target),
                call. = FALSE))
  }

  ## Note: the vulnerability curve (stem_c/stem_b/psi_crit) is derived from K_s,
  ## and TF24_cost_scale from rho, inside TF24_hyperpar (#548). Those descriptor columns are
  ## therefore intentionally absent from scenario_mapping.csv and left to the
  ## hyperpar; passing them as input traits would trip its overwrite guard.

  expected <- switch(as.character(row$Expectation),
                     "Model failure"           = "failure",
                     "Model runs successfully" = "success",
                     stop(sprintf("Unrecognised Expectation '%s' (scenario '%s').",
                                  row$Expectation, row$Scenario), call. = FALSE))

  list(traits = unlist(traits), env = env, driver = driver, expected = expected)
}

##' @return \code{scenario_control} returns the \code{Control} the gateway runs
##'   under: \code{control()} with \code{node_density_in_birth_date = TRUE}.
##'
##'   Every scenario here is TF24, and TF24 is the model the two density
##'   coordinates genuinely disagree on (#590). The transport equation's
##'   compression term is the total derivative of growth along a cohort's own
##'   trajectory, which equals \code{dg/dh} only when growth is a function of
##'   size; TF24's reserve gate (#517) breaks that, so in height coordinates the
##'   solver carries an accurate derivative of a quantity the plant never
##'   experiences. The tell is that refining the node schedule does not close the
##'   height-vs-birth-date gap the way it does for FF16 and K93, where growth
##'   *is* size-only and the two coordinates converge at ~2nd order.
##'
##'   Measured on this gateway at \code{max_patch_lifetime = 100}, the
##'   coordinate change raises R0 on every scenario, by 2.4x (S02) to 47x (S07),
##'   and moves S01 across R0 = 1. So the choice is not cosmetic, and running the
##'   TF24 gateway in the package default coordinate would score the model on a
##'   compression term that is wrong for it.
##' @rdname scenario_eval
##' @export
scenario_control <- function() {
  control(node_density_in_birth_date = TRUE)
}

##' @param config A config list from \code{scenario_to_config}.
##' @param max_patch_lifetime Patch lifetime (years) for the SCM run.
##' @param ctrl A \code{Control} object. Defaults to \code{scenario_control()},
##'   which integrates in birth date rather than height -- see there for why
##'   that is the right coordinate for TF24.
##' @param birth_rate Birth rate passed to \code{add_strategies}.
##' @return \code{build_scenario} returns a list with the \code{Parameters}
##'   (\code{p}), configured \code{Environment} (\code{env}) and \code{Control}
##'   (\code{ctrl}) ready for \code{run_scm}.
##' @rdname scenario_eval
##' @export
build_scenario <- function(config, max_patch_lifetime = 100,
                           ctrl = scenario_control(), birth_rate = 1) {
  p <- scm_base_parameters("TF24")
  p$max_patch_lifetime <- max_patch_lifetime

  tr <- config$traits
  traits <- trait_matrix(unname(tr), names(tr))
  p <- add_strategies(p, traits, hyperpar = TF24_hyperpar, birth_rate = birth_rate)

  env <- Environment("TF24")
  n_depths <- env$get_soil_number_of_depths()
  for (nm in names(config$env)) {
    env[[nm]] <- config$env[[nm]]
  }
  ## Start each layer at half saturation for a reproducible initial condition.
  env$set_soil_water_state(rep(0.428 * 0.5, n_depths))

  mean_rain <- config$driver[["rainfall_mean"]]
  amp_frac  <- config$driver[["rainfall_amp_frac"]]
  if (is.null(amp_frac)) {
    amp_frac <- 0
  }
  if (!is.null(mean_rain)) {
    if (amp_frac == 0) {
      env$extrinsic_drivers_set_constant("rainfall", mean_rain)
    } else {
      ## Sinusoidal annual rainfall; trough clamped to zero (a trough of exactly
      ## zero is the "amplitude = mean" extreme-seasonality case).
      ##
      ## Knots are placed per *year*, not per run, so the realised seasonality
      ## does not depend on max_patch_lifetime. The previous
      ## max(200, mpl * 6) gave 6 knots per annual cycle at the default
      ## mpl = 100. Measured, that is better than it looks -- cubic
      ## interpolation of a sine at 6 points/cycle is accurate to 5.4e-3 on a
      ## peak of 3.0 (0.18%), and conserves the annual total to 1e-5 % -- but it
      ## is 4th-order in the spacing and degrades with mpl, so pin the density
      ## instead. At 48/yr the error is ~1e-6, i.e. numerically exact.
      points_per_year <- 48L
      t <- seq(0, max_patch_lifetime,
               length.out = max(200L, as.integer(ceiling(max_patch_lifetime *
                                                           points_per_year))))
      y <- mean_rain * (1 + amp_frac * sin(2 * pi * t))
      y <- pmax(y, 0)
      env$extrinsic_drivers_set_variable("rainfall", t, y)
    }
  }

  list(p = p, env = env, ctrl = ctrl)
}

## Demographic persistence, as distinct from the run completing.
##
## With birth_rate = 1 (build_scenario's default) offspring_production *is* the
## net reproduction ratio R0 -- lifetime offspring per offspring introduced --
## so a strategy replaces itself only at R0 >= 1. This matters because `status`
## and `outcome` below test `total > 0`, which is a far weaker condition than it
## looks: at the blessed baseline five of the eight hydraulic scenarios return
## R0 between 2e-15 and 6e-14 -- numerically extinct -- and are nonetheless
## recorded as "persisted". That is the main reason the gateway discriminates
## so little (3/8, no expected failure failing).
##
## Reported alongside, not folded into, `status`/`outcome`: the scenario CSV's
## "Model failure" means the model *breaks numerically* (#549/#550), not that
## the strategy dies out, so the two axes are different questions and the
## blessed baseline diff is defined on the existing ones.
persists_at <- function(total, finite, threshold = 1) {
  isTRUE(finite) && isTRUE(total >= threshold)
}

##' @param p A built \code{Parameters} object (e.g. \code{build_scenario()$p}).
##' @param env An \code{Environment} object (e.g. \code{build_scenario()$env}).
##' @return \code{classify_scm_run} returns a list describing the run:
##'   \code{status} (\code{"success"}/\code{"failure"}),
##'   \code{offspring_production}, \code{persists}, \code{finite},
##'   \code{error_message}, \code{warnings} and \code{run_seconds}. A run is a
##'   success when it completes with finite, positive total offspring
##'   production; any thrown error or non-finite output is a failure. The C++
##'   layer already fails fast on non-finite state, so classification does not
##'   depend on error wording. \code{persists} is the separate, stricter
##'   ecological question — whether the strategy replaces itself, R0 >= 1 —
##'   which a run can fail while still counting as a numerical success.
##' @rdname scenario_eval
##' @export
classify_scm_run <- function(p, env, ctrl = scenario_control()) {
  ## The test suite sets options(warn = 2); make sure warnings raised during a
  ## scenario run are recorded, not escalated to errors.
  withr::local_options(warn = 1)

  warns <- character(0)
  t0 <- proc.time()[["elapsed"]]
  res <- tryCatch(
    withCallingHandlers(
      {
        scm <- run_scm(p, env = env, ctrl = ctrl)
        list(ok = TRUE, offspring_production = scm$offspring_production)
      },
      warning = function(w) {
        warns <<- c(warns, conditionMessage(w))
        invokeRestart("muffleWarning")
      }),
    error = function(e) list(ok = FALSE, error_message = conditionMessage(e)))
  elapsed <- proc.time()[["elapsed"]] - t0

  warn_str <- if (length(warns)) paste(unique(warns), collapse = " | ") else NA_character_

  if (isTRUE(res$ok)) {
    op <- res$offspring_production
    finite <- length(op) > 0L && all(is.finite(op))
    total <- if (finite) sum(op) else NA_real_
    status <- if (finite && isTRUE(total > 0)) "success" else "failure"
    ## outcome distinguishes a numerical failure (crash / non-finite) from a
    ## finite run in which the strategy simply fails to persist (extinct).
    outcome <- if (!finite) "crashed" else if (total > 0) "persisted" else "extinct"
    list(status = status, outcome = outcome, crashed = !finite,
         offspring_production = total, persists = persists_at(total, finite),
         finite = finite,
         error_message = NA_character_, warnings = warn_str, run_seconds = elapsed)
  } else {
    list(status = "failure", outcome = "crashed", crashed = TRUE,
         offspring_production = NA_real_, persists = FALSE, finite = FALSE,
         error_message = res$error_message, warnings = warn_str,
         run_seconds = elapsed)
  }
}

##' @return \code{evaluate_scenario} returns a one-row scorecard tibble.
##' @rdname scenario_eval
##' @export
evaluate_scenario <- function(row, mapping, ctrl = scenario_control(),
                              max_patch_lifetime = 100) {
  config <- scenario_to_config(row, mapping)
  built <- build_scenario(config, max_patch_lifetime = max_patch_lifetime,
                          ctrl = ctrl)
  run <- classify_scm_run(built$p, built$env, built$ctrl)

  tibble::tibble(
    scenario_id          = row$scenario_id %||% NA_character_,
    scenario             = row$Scenario,
    expected             = config$expected,
    observed             = run$status,
    match                = identical(run$status, config$expected),
    outcome              = run$outcome,
    crashed              = run$crashed,
    offspring_production = run$offspring_production,
    persists             = run$persists,
    finite               = run$finite,
    error_message        = run$error_message,
    warnings             = run$warnings,
    run_seconds          = run$run_seconds,
    config               = list(config))
}

##' @param scenarios A scenario tibble from \code{read_scenario_table}.
##' @param workers Number of parallel workers. \code{> 1} uses **fork-based**
##'   parallelism (\code{parallel::mclapply}). Forking is deliberate: it inherits
##'   the currently-loaded namespace and compiled library, so it works when
##'   \code{plant} is loaded for development via \code{pkgload::load_all} /
##'   \code{devtools::load_all}. A PSOCK / \code{future::multisession} cluster
##'   would spawn fresh R sessions that see only the *installed* package, not the
##'   dev build, and would silently run stale (or missing) code — so it is not
##'   used here. Forking is unavailable on Windows, where the run falls back to
##'   sequential regardless of \code{workers}.
##' @return \code{run_scenarios} returns a scorecard tibble (one row per
##'   scenario) with a \code{"metadata"} attribute recording provenance. A crash
##'   in one scenario never aborts the batch. Scenario runs are deterministic
##'   (no RNG), so parallel and sequential runs produce identical scorecards.
##' @param cache Optional path to an \code{.rds} cache. When supplied, each
##'   scenario's result is keyed by a content hash of its resolved config,
##'   \code{max_patch_lifetime} and the model fingerprint
##'   (\code{\link{scenario_model_fingerprint}}); a scenario is rerun only when
##'   that key changes. The key covers every input that affects a result: the
##'   resolved config, \code{max_patch_lifetime}, all \code{ctrl} settings, and
##'   the model fingerprint. So a model recompile / R-source edit reruns
##'   everything (the fingerprint moves), a \code{ctrl} change reruns everything,
##'   and editing one mapping cell reruns only the scenarios it touches. The
##'   cache deliberately errs toward rerunning.
##' @rdname scenario_eval
##' @export
run_scenarios <- function(scenarios = read_scenario_table(),
                          mapping = read_scenario_mapping(),
                          ctrl = scenario_control(), max_patch_lifetime = 100,
                          workers = 1L, cache = NULL) {
  eval_error_row <- function(row, msg) tibble::tibble(
    scenario_id = row$scenario_id %||% NA_character_,
    scenario = row$Scenario, expected = NA_character_,
    observed = "error", match = NA, outcome = "error", crashed = TRUE,
    offspring_production = NA_real_, persists = FALSE, finite = FALSE,
    error_message = msg, warnings = NA_character_,
    run_seconds = NA_real_, config = list(NULL))

  eval_one <- function(i) {
    row <- scenarios[i, , drop = FALSE]
    tryCatch(
      evaluate_scenario(row, mapping, ctrl = ctrl,
                        max_patch_lifetime = max_patch_lifetime),
      error = function(e) eval_error_row(row, conditionMessage(e)))
  }

  n <- nrow(scenarios)
  idx <- seq_len(n)

  ## Content-addressed key per scenario: every input that affects the result —
  ## resolved config, lifetime, the Control settings, and the model fingerprint.
  keys <- NULL
  cached <- NULL
  if (!is.null(cache)) {
    fp <- scenario_model_fingerprint()
    ctrl_vals <- control_values(ctrl)
    keys <- vapply(idx, function(i) {
      cfg <- tryCatch(scenario_to_config(scenarios[i, , drop = FALSE], mapping),
                      error = function(e) NULL)
      rlang::hash(list(config = cfg, mpl = max_patch_lifetime,
                       ctrl = ctrl_vals, fingerprint = fp))
    }, character(1))
    if (file.exists(cache)) {
      cached <- tryCatch(readRDS(cache), error = function(e) NULL)
    }
  }

  ## Which scenarios can be reused from the cache (matching key)?
  reuse_row <- vector("list", n)
  cached_keys <- if (!is.null(cached)) attr(cached, "keys") else NULL
  if (!is.null(keys) && !is.null(cached_keys)) {
    hit <- match(keys, cached_keys)
    for (i in idx[!is.na(hit)]) {
      row <- cached[hit[i], , drop = FALSE]
      ## Refresh identity columns in case scenario order/names changed.
      row$scenario_id <- scenarios$scenario_id[i] %||% row$scenario_id
      row$scenario <- scenarios$Scenario[i]
      reuse_row[[i]] <- row
    }
  }
  run_idx <- idx[vapply(reuse_row, is.null, logical(1))]
  if (!is.null(cache)) {
    message(sprintf("Scenario cache: %d reused, %d to run.",
                    n - length(run_idx), length(run_idx)))
  }

  use_fork <- workers > 1L && .Platform$OS.type != "windows"
  if (length(run_idx) && use_fork) {
    ## mc.preschedule = FALSE load-balances the uneven per-scenario runtimes.
    fresh <- parallel::mclapply(run_idx, eval_one, mc.cores = workers,
                                mc.preschedule = FALSE)
    ## A worker that died (rather than returning a row) surfaces as a
    ## try-error; turn it into a visible error row rather than dropping it.
    fresh <- Map(function(res, i) {
      if (inherits(res, "try-error") || !tibble::is_tibble(res)) {
        eval_error_row(scenarios[i, , drop = FALSE],
                       paste("worker failed:", as.character(res)))
      } else {
        res
      }
    }, fresh, run_idx)
  } else {
    fresh <- lapply(run_idx, eval_one)
  }

  rows <- reuse_row
  rows[run_idx] <- fresh
  scorecard <- dplyr::bind_rows(rows)

  meta <- scenario_run_metadata()
  meta$workers <- if (use_fork) workers else 1L
  ## Record the density coordinate and the patch lifetime. Both change the
  ## numbers -- the coordinate by up to 47x, enough to flip a persistence
  ## verdict -- so a stored scorecard that does not say which produced it cannot
  ## honestly be compared against another, and the blessed baseline is exactly
  ## such a stored scorecard.
  meta$node_density_in_birth_date <- ctrl$node_density_in_birth_date
  meta$max_patch_lifetime <- max_patch_lifetime
  attr(scorecard, "metadata") <- meta
  attr(scorecard, "max_patch_lifetime") <- max_patch_lifetime
  if (!is.null(keys)) {
    attr(scorecard, "keys") <- keys
    if (!is.null(cache)) {
      saveRDS(scorecard, cache)
    }
  }
  scorecard
}

##' @return \code{scenario_model_fingerprint} returns a hash string that changes
##'   whenever the model that produces scenario results could change: the
##'   package version, the compiled shared library, and all package R source
##'   plus the scenario data files. It is intentionally broad — the cache should
##'   rerun rather than risk a stale result.
##' @rdname scenario_eval
##' @export
scenario_model_fingerprint <- function() {
  ## Compiled C++ (covers every strategy/environment/solver change).
  dll <- tryCatch(getLoadedDLLs()[["plant"]][["path"]],
                  error = function(e) NA_character_)
  dll_hash <- if (!is.na(dll) && file.exists(dll))
    unname(tools::md5sum(dll)) else NA_character_

  ## Package R source (covers hyperpar / run_scm / engine changes) and the
  ## scenario data files. Best-effort: only hashes what is findable on disk
  ## (present in a load_all dev tree; absent for an installed package).
  src <- c(list.files("R", pattern = "\\.[Rr]$", full.names = TRUE),
           scenario_file("model_scenarios_hydraulic.csv"),
           scenario_file("scenario_mapping.csv"))
  src <- src[file.exists(src)]
  src_hashes <- if (length(src)) unname(tools::md5sum(src)) else character(0)

  rlang::hash(list(
    package_version = as.character(utils::packageVersion("plant")),
    dll = dll_hash,
    src = sort(src_hashes)))
}

## Stable snapshot of a Control object's data fields, for cache keying. Reads
## every (non-function) field by name so any control change alters the key.
control_values <- function(ctrl) {
  nm <- sort(names(ctrl))
  vals <- lapply(nm, function(n) {
    v <- tryCatch(ctrl[[n]], error = function(e) NULL)
    if (is.function(v)) NULL else v
  })
  names(vals) <- nm
  vals[!vapply(vals, is.null, logical(1))]
}

##' @return \code{scenario_run_metadata} returns a list of provenance fields
##'   (git commit / branch / dirty flag, package version, R version, platform,
##'   timestamp) so a scorecard records the exact build it came from.
##' @rdname scenario_eval
##' @export
scenario_run_metadata <- function() {
  git <- function(args) {
    tryCatch(
      {
        out <- suppressWarnings(system2("git", args, stdout = TRUE, stderr = FALSE))
        if (length(out)) trimws(paste(out, collapse = "\n")) else NA_character_
      },
      error = function(e) NA_character_)
  }
  dirty <- tryCatch(
    length(suppressWarnings(system2("git", c("status", "--porcelain"),
                                    stdout = TRUE, stderr = FALSE))) > 0L,
    error = function(e) NA)
  list(
    git_commit      = git(c("rev-parse", "HEAD")),
    git_branch      = git(c("rev-parse", "--abbrev-ref", "HEAD")),
    git_dirty       = dirty,
    package_version = as.character(utils::packageVersion("plant")),
    r_version       = R.version.string,
    platform        = R.version$platform,
    timestamp       = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
}

##' @param scorecard A scorecard tibble from \code{run_scenarios}.
##' @return \code{scenario_summary} returns a one-row tibble reporting **two
##'   separate axes** (#572), because the gateway conflated them and so returned
##'   no signal:
##'
##'   \describe{
##'     \item{Numerical viability — did the model run?}{\code{n},
##'       \code{n_ran}, \code{n_crashed}, \code{viability_rate}. This is what the
##'       scenario CSV's "Model failure" means (#549, #550) and what the
##'       hydraulic/NSC work targets.}
##'     \item{Ecological persistence — does the strategy replace itself?}{
##'       \code{n_persists}, \code{persistence_rate}, judged at R0 >= 1 (see
##'       \code{persists_at}). Now that the crashes are fixed this is the axis the
##'       gateway is actually useful for.}
##'   }
##'
##'   \code{n_match} / \code{match_rate} and the expected-failure vs
##'   expected-success breakdown are still reported, but as **agreement with the
##'   CSV's crash predictions**, not as a headline quality score: with the crashes
##'   fixed, the match rate mostly measures how well those predictions have aged.
##'   Read \code{n_persists} against \code{n} — a scorecard where every run
##'   "succeeds" but almost none persists is reporting that the model no longer
##'   crashes, not that the strategies live.
##' @rdname scenario_eval
##' @export
scenario_summary <- function(scorecard) {
  exp_fail <- scorecard$expected == "failure"
  exp_succ <- scorecard$expected == "success"
  n <- nrow(scorecard)
  ## Scorecards recorded before `persists` / `crashed` existed -- including the
  ## blessed baseline in tests/testthat/test_data/ -- have no such column, and
  ## must still summarise rather than error.
  persists <- if ("persists" %in% names(scorecard)) scorecard$persists else NA
  ## Numerical viability. Prefer the recorded `crashed` flag; fall back to
  ## `observed == "failure"` only for scorecards predating that column, where the
  ## two coincide because `status` was the only axis.
  crashed <- if ("crashed" %in% names(scorecard)) {
    scorecard$crashed
  } else {
    scorecard$observed != "success"
  }
  n_crashed <- sum(crashed, na.rm = TRUE)
  n_persists <- sum(persists, na.rm = TRUE)
  tibble::tibble(
    n                      = n,
    ## --- axis 1: numerical viability (did it run?) ---
    n_ran                  = n - n_crashed,
    n_crashed              = n_crashed,
    viability_rate         = if (n > 0) (n - n_crashed) / n else NA_real_,
    ## --- axis 2: ecological persistence (does the strategy live?) ---
    n_persists             = n_persists,
    persistence_rate       = if (n > 0) n_persists / n else NA_real_,
    ## --- agreement with the CSV's (crash-era) expectations ---
    n_match                = sum(scorecard$match, na.rm = TRUE),
    match_rate             = mean(scorecard$match, na.rm = TRUE),
    n_expected_fail        = sum(exp_fail, na.rm = TRUE),
    n_expected_fail_met    = sum(exp_fail & scorecard$match, na.rm = TRUE),
    n_expected_success     = sum(exp_succ, na.rm = TRUE),
    n_expected_success_met = sum(exp_succ & scorecard$match, na.rm = TRUE))
}

##' @param output_file Output HTML path for the rendered report.
##' @param input_file The report template (\code{.Rmd}).
##' @param overwrite Overwrite an existing \code{output_file}?
##' @param quiet Passed to \code{rmarkdown::render}.
##' @return \code{scenario_generate_report} renders the scorecard report and
##'   returns the output path (invisibly).
##' @rdname scenario_eval
##' @export
scenario_generate_report <- function(scorecard,
                                     output_file = "scenario_scorecard.html",
                                     input_file = system.file(
                                       "reports", "scenario_scorecard.Rmd",
                                       package = "plant"),
                                     overwrite = FALSE, quiet = TRUE) {
  if (!nzchar(input_file) || !file.exists(input_file)) {
    stop("Could not find the scorecard report template.", call. = FALSE)
  }
  if (file.exists(output_file) && !overwrite) {
    stop(sprintf("'%s' exists; pass overwrite = TRUE to replace it.", output_file),
         call. = FALSE)
  }
  output_file <- normalizePath(output_file, mustWork = FALSE)
  rmarkdown::render(input_file, output_file = output_file,
                    params = list(scorecard = scorecard), quiet = quiet,
                    envir = new.env(parent = globalenv()))
  invisible(output_file)
}

## Local null-coalescing helper (avoids taking a hard rlang dependency here).
`%||%` <- function(a, b) if (is.null(a)) b else a
