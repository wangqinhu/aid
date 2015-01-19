#' Analysis of inoculation data
#'
#' Analysis and illustration of inoculation data.
#' 
#' @param file a text file contains inoculation data.
#' @param type grade, lesion or biomass (qPCR), must be specified
#' @param alternative a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less". You can specify just the initial letter.
#' @return statistical analysis and illustration for inoculation data
#' @seealso \code{\link{grade.test}}, \code{\link{lesion.test}}, \code{\link{biomass.test}}, \code{\link{grade.barplot}}, \code{\link{lesion.barplot}} and \code{\link{biomass.barplot}}.
#' @export
#' @examples
#' library(aid)
#' # grade data
#' demo1 <- system.file("extdata", "demo1.tsv", package="aid")
#' aid(demo1, type = "grade", alternative="less")
#' # lesion data
#' demo2 <- system.file("extdata", "demo2.tsv", package="aid")
#' aid(demo2, type = "lesion", alternative="two.sided")
#' demo3 <- system.file("extdata", "demo3.tsv", package="aid")
#' aid(demo3, type = "biomass", alternative="two.sided")
aid <- function (file, type, alternative="two.sided") {
  ino <- read.table(file, header = TRUE, check.names=FALSE)
  if (type == "grade") {
    grade.barplot(ino, alternative)
  } else if (type == "lesion") {
    lesion.barplot(ino, alternative)
  } else if (type == "biomass") {
    biomass.barplot(ino, alternative)
  } else {
    stop("Unknown inoculation data type!\n")
  }
}