#' Statistical test for grade data
#'
#' Convert grade data to real observations and perform t test.
#' 
#' @param ino a list contains inoculation grade data.
#' @param alternative a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less". You can specify just the initial letter.
#' @return p-value symbols indicating the significances between comparisons.
#' @seealso \code{\link{aid}}, \code{\link{grade.barplot}}, \code{\link{lesion.test}} and \code{\link{biomass.test}}.
#' @export
#' @examples
#' demo <- system.file("extdata", "demo1.tsv", package="aid")
#' dat <- read.table(demo, header = TRUE, check.names=FALSE)
#' grade.test(dat, alternative = "t")
grade.test <- function(ino, alternative) {
  
  # number of individuals and grades
  ni <- dim(ino)[1]
  ng <- dim(ino)[2]
  
  # initialization of p-value and p-value signs 
  ino.p <- rep(1, ni)
  ino.ps <- rep(" ", ni)
  
  # for the control individual
  a <- NULL
  for (j in 1:ng) {
    if (ino[1,j] > 0) {
      a[(length(a) + 1):(length(a) + ino[1,j])] <- rep(j, ino[1,j])
    }
  }
  
  for (i in 2:ni) {   
    # for the treatment individual
    b <- NULL
    for (j in 1:ng) {
      if (ino[i,j] > 0) {
        b[(length(b)+1):(length(b) + ino[i,j])] <- rep(j, ino[i,j])
      }
    }
    
    # perform t test
    t <- t.test(a, b, alternative)
    
    # assign p-value symbols
    ino.p[i] <- t$p.value
    ino.ps[i] <- sym.pval(ino.p[i])
    
    i <- i + 1;
    
  }
  
  return(ino.ps)
  
}


#' Statistical test for lesion data
#'
#' Perform t test on lesion data.
#'
#' @param ino a list contains inoculation lesion data.
#' @param alternative a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less". You can specify just the initial letter.
#' @return p-value symbols indicating the significances between comparisons.
#' @seealso \code{\link{aid}}, \code{\link{lesion.barplot}}, \code{\link{biomass.test}} and \code{\link{grade.test}}.
#' @export
#' @examples
#' demo <- system.file("extdata", "demo2.tsv", package="aid")
#' dat <- read.table(demo, header = TRUE, check.names=FALSE)
#' lesion.test(dat, alternative = "t")
lesion.test <- function(ino, alternative) {
  
  # number of individuals
  ni <- length(ino)
  
  # initialization of p-value and p-value signs
  ino.p <- rep(1, ni)
  ino.ps <- rep(" ", ni)
  
  for (i in 2:ni) {
    t <- t.test(ino[,1], ino[,i], alternative)
    ino.p[i] <- t$p.value
    ino.ps[i] <- sym.pval(ino.p[i])
    i < i + 1
  }
  
  return(ino.ps)
  
}


#' Statistical test for biomass data
#'
#' Perform t test on biomass (qPCR) data.
#'
#' @param ino a list contains inoculation biomass data.
#' @param alternative a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less". You can specify just the initial letter.
#' @return p-value symbols indicating the significances between comparisons.
#' @seealso \code{\link{aid}}, \code{\link{biomass.barplot}}, \code{\link{lesion.test}} and \code{\link{grade.test}}.
#' @export
#' @examples
#' bio<-c(0.82, 3.14, 0.88, 3.21, 0.85, 3.20)
#' dim(bio)<-c(2,3)
#' biomass.test(bio, alternative = "t")
biomass.test <- function(bio, alternative) {
  
  # number of individuals
  ni <- dim(bio)[1]
  
  # initialization of p-value and p-value signs
  ino.p <- rep(1, ni)
  ino.ps <- rep(" ", ni)
  
  for (i in 2:ni) {
    t <- t.test(bio[1,], bio[i,], alternative)
    ino.p[i] <- t$p.value
    ino.ps[i] <- sym.pval(ino.p[i])
    i < i + 1
  }
  
  return(ino.ps)
  
}