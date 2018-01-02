table_plant_parameters <- function(filename, dest){
  tab <- read.csv(filename, stringsAsFactors=FALSE,
                  check.names=FALSE, comment.char="#")

  to_bold <- tab$Symbol == ""
  tab$Description[to_bold] <- sprintf("\\textbf{%s}", tab$Description[to_bold])

  to_math <- tab$Symbol != ""
  tab$Symbol[to_math] <- sprintf("$%s$", tab$Symbol[to_math])

  to_math <- grepl("[_^]", tab$Unit)
  tab$Unit[to_math] <-
    sprintf("$%s$",
            gsub("(mol|kg|m|yr|d)", "\\\\mathrm{\\1}", tab$Unit[to_math]))
  tab$Unit[to_math] <-
            gsub(" ", "\\,", tab$Unit[to_math], fixed =TRUE)

  s <- c(FF16_Strategy(),
         lapply(as.list(args(make_FF16_hyperpar)), eval))
  tab$Value <- NA_real_

  # Now extract actual values from plant
  to_code <- tab[["Code"]] != "" & tab[["Code"]] %in% names(s)
  code <- tab[["Code"]][to_code]

  val <- s[code]
  err <- sapply(val, is.null)
  if (any(err)) {
    message("Missing parameters: ", paste(code[err], collapse=", "))
  }

  oo <- options(scipen=999)
  on.exit(options(oo))
  tab$Value[to_code] <- unlist(lapply(val, prettyNum))

  # Display names in latex
  to_code <- tab[["Code"]] != ""

  tab[["Code"]][to_code] <- sprintf("\\texttt{%s}",
                                      gsub("_", "\\\\_", tab[["Code"]][to_code]))

  x <- xtable(tab,
              hline.after=c(1),
              align='lp{7cm}llll')
  y <- print(x, sanitize.text.function=as.character,
             include.rownames=FALSE, floating=FALSE,
             print.results=FALSE)
  cat(y, file=dest)
}
