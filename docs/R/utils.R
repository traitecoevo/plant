empty_box <- function() {
  par(oma=c(0,0,0,0), mar=rep(0.1,4))
  plot(1,1, ann=FALSE, axes=FALSE, type='n')
  box()
}

label <- function(text, px=0.03, py=NULL, ..., adj=c(0, 1)) {
  if (is.null(py)) {
    fin <- par("fin")
    r <- fin[[1]] / fin[[2]]
    if (r > 1) { # x is longer.
      py <- 1 - px
      px <- (1 - py) / r
    } else {
      py <- 1 - px * r
    }
  }
  usr <- par("usr")
  x <- usr[1] + px*(usr[2] - usr[1])
  y <- usr[3] + py*(usr[4] - usr[3])

  ## NOTE: base 10 log:
  if (par("xlog")) {
    x <- 10^x
  }
  if (par("ylog")) {
    y <- 10^y
  }

  text(x, y, text, adj=adj, ...)
}

label_panel <- function(i, ...) {
  label(sprintf("%s)", letters[i]), ...)
}

make_transparent <- function(col, opacity=.5) {
  alpha <- opacity
  if (length(alpha) > 1 && any(is.na(alpha))) {
    n <- max(length(col), length(alpha))
    alpha <- rep(alpha, length.out=n)
    col <- rep(col, length.out=n)
    ok <- !is.na(alpha)
    ret <- rep(NA, length(col))
    ret[ok] <- make_transparent(col[ok], alpha[ok])
    ret
  } else {
    tmp <- col2rgb(col)/255
    rgb(tmp[1,], tmp[2,], tmp[3,], alpha=alpha)
  }
}

mix_colours <- function(col1, col2, p) {
  m1 <- drop(col2rgb(col1))
  m2 <- drop(col2rgb(col2))
  m3 <- outer(m1, as.numeric(p)) + outer(m2, as.numeric(1 - p))
  ret <- rep(NA_character_, length(p))
  i <- colSums(is.na(m3)) == 0L
  ret[i] <- rgb(m3[1, i], m3[2, i], m3[3, i], maxColorValue=255)
  ret
}

## no-extrapolation version of splinefun_loglog
splinefun_loglog2 <- function(x, y, xx){
  ret <- splinefun_loglog(x, y)(xx)
  ii <- xx <= max(x, na.rm=TRUE) & xx >= min(x, na.rm=TRUE)
  ret[!ii] <- 0
  ret
}

combine_md <- function(..., output) {
  find_metadata <- function(x) {
    i <- grep("^---+\\s*$", x)
    if (length(i) < 2L) {
      stop("Metadata block not found")
    }
    if (i[[1]] != 1L) {
      stop("Metadata block not the first line")
    }
    i
  }
  demote <- function(x) {
    i <- which(grepl("^#+", x))
    i <- i[i == 1 | grepl("^\\s*$", x[i - 1])]
    x[i] <- paste0("#", x[i])
    x
  }

  files <- c(...)
  dat <- lapply(files, readLines)
  titles <- character(length(dat))

  for (i in seq_along(dat)) {
    x <- dat[[i]]
    j <- find_metadata(x)
    t <- yaml::yaml.load(x[(j[[1]] + 1L):(j[[2]] - 1L)])$title
    t <- sub("plant:[^:]*: ", "", t)
    t <- gsub("(^_|_$)", "", t)
    titles[[i]] <- t
    dat[[i]] <- c("", paste("#", t), demote(x[-(j[[1]]:j[[2]])]))
  }

  header <-
    c("---",
      'title: "Worked examples showing how to interact with plant from R"',
      "---")
  ret <- c(header, unlist(dat))
  writeLines(ret, output)
}

## I'm 99% sure we've done this before.
md_to_pdf <- function(filename) {
  output <- sub("\\.md", ".pdf", basename(filename))
  template <- c("--template", "template.tex")
  opts <- "--toc"
  owd <- setwd(dirname(filename))
  on.exit(setwd(owd))
  call_system(Sys_which("pandoc"),
              c(basename(filename), "-o", output, template, opts))
}
