
# Causes pkgdown to use bookdown::html_document2 output, as this enables referencing
# to figures, equations and tables
# From https://gist.github.com/hturner/3152081e223ade0bb212bcef19f183bf
# code by Andrzej K. Oleś @aoles
# https://github.com/hadley/pkgdown/issues/323#issuecomment-341907743

library(pkgdown)
assignInNamespace(
    "build_rmarkdown_format",
    ns = "pkgdown",
    value = function(pkg = ".",
                     depth = 1L,
                     data = list(),
                     toc = TRUE) {
        # Render vignette template to temporary file
        path <- tempfile(fileext = ".html")
        suppressMessages(pkgdown::render_page(pkg, "vignette", data, path, depth = depth))
        
        list(
            path = path,
            format = bookdown::html_document2(
                toc = toc,
                toc_depth = 2,
                self_contained = FALSE,
                theme = NULL,
                template = path,
                number_sections = FALSE
            )
        )
    }
)
