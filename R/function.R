#' Assign p-value symbols
#'
#' Convert p-value to symbols
#'
#' @param pval pvalue, between 0-1.
#' @seealso \code{\link{aid}}, \code{\link{grade.test}}, \code{\link{lesion.test}} and \code{\link{biomass.test}}.
#' @export
#' @examples
#' pval <- c(0.15, 0.10, 0.05, 0.02, 0.01, 0.0001)
#' n <- length(pval)
#' sym  <- rep(" ", n)
#' for (i in 1:n)
#'   sym[i] <- sym.pval(pval[i])
#' cat(sym)
sym.pval <- function(pval) {
  
  if (pval < 0 || pval > 1)
    stop("p-value must be between 0 and 1 !\n")
  
  # threshold
  th = c(0.001, 0.01, 0.05, 0.1)
  # symbols
  sy = c("***", "**", "*", "+", " ")
  
  if (pval <= th[1]) {
    ps <- sy[1]
  } else if (pval <= th[2]) {
    ps <- sy[2]
  } else if (pval <= th[3]) {
    ps <- sy[3]
  } else if (pval <= th[4]) {
    ps <- sy[4]
  } else {
    ps <- sy[5]
  }
  
  return(ps)
  
}


# Disease severity index
#'
#' Caculate disease severity
#' 
#' @param ino a list contains inoculation grade data.
#' @return Disease severity index for individuals.
#' @seealso \code{\link{aid}}, \code{\link{grade.barplot}} and \code{\link{grade.test}}.
#' @export
#' @examples
#' demo <- system.file("extdata", "demo1.tsv", package="aid")
#' dat <- read.table(demo, header = TRUE, check.names=FALSE)
#' dsi(dat)
dsi <- function(ino) {
  
  # number of individuals and grades
  ni <- dim(ino)[1]
  ng <- dim(ino)[2]
  
  # grades
  gr <- as.numeric(colnames(ino))
  
  # initialization of dsi
  ino.dsi <- rep(0, ni)
    
  for (i in 1:ni) {
    for (j in 1:ng) {
      ino.dsi[i] <- formatC( sum( gr[j] * ino[i, j] / sum(ino[i,]) ), digits=1 )
    }
  }
  
  return(ino.dsi)
  
}