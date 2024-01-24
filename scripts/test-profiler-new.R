clean_search <- function() {
  defaults <- c(".GlobalEnv",
                paste0("package:", getOption("defaultPackages")),
                "Autoloads",
                "package:base")
  currentList <- search()
  deletes <- setdiff(currentList, defaults)
  for (entry in deletes)
    detach(entry, character.only = TRUE)
}

# function to knit directory
knit_dir <- function(dir) {
  # create the full path
  full_dir <- normalizePath(dir)
  
  # list the Rmd files with full names
  files <-
    list.files(full_dir, pattern = '*.Rmd$', full.names = TRUE)
  
  # render the files using all output formats in the YAML
  purrr::map(files, function(file) {
    # clean search list to avoid conflicts
    clean_search()
    
    # render the file using all
    # output formats in a new env
    rmarkdown::render(file,
                      output_format = "all",
                      envir = new.env())
  })
}


test_that("offspring arrival", {
  out_file <- tempfile("jointprof", fileext = ".pb.gz")
  jointprof::start_profiler(out_file)
  
    p0 <- scm_base_parameters("FF16w")
    p0$max_patch_lifetime = 0.433
    
    p1 <- expand_parameters(trait_matrix(c(2), "p_50"), p0, mutant=FALSE, birth_rate_list = c(1))
    
    env <- make_environment("FF16w", soil_initial_state = rep(0.2, 1), rainfall = 1.5)
    
    ctrl <- scm_base_control()
    
    for (i in 1:1000) {
    out <- run_scm(p1, env, ctrl) 
  }
  
  profile_data <- jointprof::stop_profiler()
  pprof_file <- tempfile("jointprof", fileext = ".pb.gz")
  profile::write_pprof(profile_data, pprof_file)
  dir.create("jointprof_fig", recursive = TRUE, showWarnings = FALSE)
  svg_file <- "jointprof_fig/minimal.svg"
  
  system2(
    jointprof::find_pprof(),
    c(
      "-svg",
      "-nodefraction 0.01",
      "-output",
      shQuote(svg_file),
      shQuote(pprof_file)
    )
  )
  
  png_file <- "jointprof_fig/minimal.png"
  rsvg::rsvg_png(svg_file, png_file, width = 5000)
  expect_equal(out$patch$environment$soil$states,
               c(22.04, -8.71, -17.59),
               tolerance=1e-02)
})
