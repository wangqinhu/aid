#' The inoculation colors vector
#'
#' Create a color vector.
#' 
#' @param n number of colors, usually equal to (or greater that) the number of desease grades.
#' @param alpha the alpha transparency, a number in [0,1]
#' @return a vector of contiguous colors.
#' @seealso \code{\link{aid}} and \code{\link{grade.barplot}}.
#' @export
#' @examples
#' ino.colors(10)
#' plot(1:5, pch = 16, col = ino.colors(5), cex=3)
ino.colors <- function (n, alpha = 1) {
  if ((n <- as.integer(n[1L])) > 0) {
    k <- n%/%2
    h <- c(4/12, 2/12, 0/12)
    s <- c(1, 0.5, 1)
    v <- c(1, 0.5, 1)
    c(hsv(h = seq.int(h[1L], h[2L], length.out = k),
          s = seq.int(s[1L], s[2L], length.out = k),
          v = seq.int(v[1L], v[2L], length.out = k),
          alpha = alpha),
      hsv(h = seq.int(h[2L], h[3L], length.out = n - k + 1)[-1L],
          s = seq.int(s[2L],  s[3L], length.out = n - k + 1)[-1L],
          v = seq.int(v[2L], v[3L], length.out = n - k + 1)[-1L],
          alpha = alpha))
  }
  else character()
}

#' Statistical test for grade data
#'
#' Convert grade data to real observations and perform t test.
#' 
#' @param ino a list contains inoculation grade data.
#' @return p-value signs indicating the significances between comparisons.
#' @seealso \code{\link{aid}}, and \code{\link{grade.barplot}}.
#' @export
#' @examples
#' demo <- system.file("extdata", "demo1.tsv", package="aid")
#' dat <- read.table(demo, header = TRUE, check.names=FALSE)
#' grade.test(dat)
grade.test <- function(ino) {

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
    t <- t.test(a, b)
    
    # assign p-value symbols
    ino.p[i] <- t$p.value
    ino.ps[i] <- sym.pval(ino.p[i])

    i <- i + 1;

  }

  return(ino.ps)

}

#' Plot inoculation data
#'
#'
#'
#'
#' Create a barplot for grade inoculation data
#' 
#' @param  ino a list contains inoculation grade data.
#' @seealso \code{\link{aid}} and \code{\link{ino.colors}}.
#' @export
#' @examples
#' demo <- system.file("extdata", "demo1.tsv", package="aid")
#' dat <- read.table(demo, header = TRUE, check.names=FALSE)
#' grade.barplot(dat)
grade.barplot <- function(ino) {

  # number of grades
  ng <- length(ino)

  # perform t test
  ps <- grade.test(ino)

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

# assign p-value symbols
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

#Error bar fuction
error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}


file <- system.file("extdata", "demo2.tsv", package="aid")
ino <- read.table(file, header = TRUE, check.names=FALSE)

lesion.test <- function(ino) {

  # number of individuals
  ni <- length(ino)

  # initialization of p-value and p-value signs
  ino.p <- rep(1, ni)
  ino.ps <- rep(" ", ni)

  for (i in 2:ni) {
    t <- t.test(ino[,1], ino[,i])
    ino.p[i] <- t$p.value
    ino.ps[i] <- sym.pval(ino.p[i])
    i < i + 1
  }

  return(ino.ps)

}

# Lesion size
lesion.barplot <- function(ino) {

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
  ps <- lesion.test(ino)
  text(barx, h + 0.5, ps)

}


#' Analysis of inoculation data
#'
#' Analysis and illustration of inoculation data.
#' 
#' @param file a text file contains inoculation data.
#' @return statistical analysis and illustration for inoculation data
#' @seealso \code{\link{grade.test}} and \code{\link{grade.barplot}}.
#' @export
#' @examples
#' library(aid)
#' demo <- system.file("extdata", "demo1.tsv", package="aid")
#' aid(demo)
aid <- function (file, type) {
  ino <- read.table(file, header = TRUE, check.names=FALSE)
  if (type == "grade") {
    grade.barplot(ino)
  } else if (type == "lesion") {
    lesion.barplot(ino)
  } else {
    stop("Unknown inoculation data type!\n")
  }
}
