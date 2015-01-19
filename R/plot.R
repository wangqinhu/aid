#' Plot error bar
#'
#' Adding error bar to a barplot
#'
#' @param x position on x-axis.
#' @param y position on y-axis.
#' @param upper length of upper error bar
#' @seealso \code{\link{aid}}, \code{\link{lesion.barplot}} and \code{\link{biomass.barplot}}.
#' @export
#' @examples
#' barx <- barplot(1:10, ylim=c(0,12))
#' error.bar(barx, 1:10, rep(1,10))
error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}


#' Plot grade data
#'
#' Create a barplot for grade inoculation data
#' 
#' @param ino a list contains inoculation grade data.
#' @param alternative a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less". You can specify just the initial letter.
#' @seealso \code{\link{aid}}, \code{\link{ino.colors}}, \code{\link{grade.test}}, \code{\link{lesion.barplot}} and \code{\link{biomass.barplot}}..
#' @export
#' @examples
#' demo <- system.file("extdata", "demo1.tsv", package="aid")
#' dat <- read.table(demo, header = TRUE, check.names=FALSE)
#' grade.barplot(dat, alternative)
grade.barplot <- function(ino, alternative) {
  
  # number of grades
  ng <- length(ino)
  
  # perform t test
  ps <- grade.test(ino, alternative)
  
  # plot
  barx <- barplot(prop.table(t(as.matrix(ino)),2),
                  col=ino.colors(ng),
                  ylim = c(0,1.25),
                  #hor
                  axes = FALSE)
  
  # add total number of leaves
  text(barx, 1.05, rowSums(ino))
  # add significant signs
  text(barx, 1.10, ps)
  # add legend
  legend("top", legend = colnames(ino),
         horiz = TRUE,
         fill=ino.colors(ng), box.col = "white")
  
}


#' Plot lesion data
#'
#' Create a barplot for lesion inoculation data
#'
#' @param ino a list contains inoculation lesion data.
#' @param alternative a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less". You can specify just the initial letter.
#' @seealso \code{\link{aid}}, \code{\link{lesion.test}}, \code{\link{biomass.barplot}} and \code{\link{grade.barplot}}.
#' @export
#' @examples
#' demo <- system.file("extdata", "demo2.tsv", package="aid")
#' dat <- read.table(demo, header = TRUE, check.names=FALSE)
#' lesion.barplot(dat, alternative)
lesion.barplot <- function(ino, alternative) {
  
  ino.mean <- colMeans(ino)
  # number of individuals
  ni <- length(ino)
  ino.sd <- NULL
  for (i in 1: ni) {
    ino.sd[i] <- sd(ino[,i])
  }
  h <-  ino.mean + ino.sd
  ymax <- max(h) + 1
  barx <- barplot(ino.mean, col=1,
                  ylim=c(0, ymax),
                  ylab="Lesion size")
  error.bar(barx, ino.mean, ino.sd)
  ps <- lesion.test(ino, alternative)
  text(barx, h + 0.5, ps)
  
}


#' Plot biomass data
#'
#' Create a barplot for biomass (qPCR) inoculation data
#'
#' @param  ino a list contains inoculation biomass data.
#' @param alternative a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less". You can specify just the initial letter.
#' @seealso \code{\link{aid}}, \code{\link{biomass.test}}, \code{\link{lesion.barplot}} and \code{\link{grade.barplot}}.
#' @export
#' @examples
#' demo <- system.file("extdata", "demo3.tsv", package="aid")
#' dat <- read.table(demo, header = TRUE, check.names=FALSE)
#' biomass.barplot(dat, alternative)
biomass.barplot <- function(ino, alternative) {
  
  # qPCR ct
  ct <- ino
  # Number of sample
  num_sam <- dim(ino)[1]
  # Number of repeat (biological or technical)
  num_rep <- dim(ino)[2]/2
  # Line of control
  lctrl <- 1
  
  # calculate relative biomass by ddct method for qPCR
  bio<-rep(NA, num_sam * num_rep)
  dim(bio)<-c(num_sam, num_rep)
  
  # ctr_ref
  ref_calibrator<-mean(as.numeric(ct[lctrl,1:num_rep]))
  calibrator<-mean(as.numeric(ct[lctrl,(num_rep+1):(2*num_rep)]-ref_calibrator))
  # bio <- 2^ddct
  for (i in 1:num_sam) {
    ref<-mean(as.numeric(ct[i,1:num_rep]))
    # dCt
    dct<-ct[i,(num_rep+1):(2*num_rep)]-ref
    # ddCt
    ddct<-dct-calibrator
    # fold
    bio[i,1:num_rep]<-2^-ddct
  }
  
  # fold
  fold<-t(bio)
  
  fold.means=rep(NA, num_sam)
  fold.sd=rep(NA, num_sam)
  
  for (i in 1:num_sam) {
    fold.means[i]<-mean(fold[,i])
    fold.sd[i]<-sd(fold[,i])
  }
  
  h <- fold.means + fold.sd
  ymax <- max(h) + 1
  barx <- barplot(fold.means, col=1, ylim=c(0,ymax+1),
                  names.arg=row.names(ct),
                  xlab="Individuals", ylab="Pathogen/Host ratio")
  error.bar(barx, fold.means, fold.sd)
  
  ps <- biomass.test(bio, alternative)
  text(barx, h + 0.5, ps)
  
}