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