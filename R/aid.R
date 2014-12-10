#' The inoculation colors vector
#'
#' Create a color vector.
#' 
#' @param number of colors, usually equal to (or greater that) the number of desease grades.
#' @return a vector of contiguous colors.
#' @seealso \code{\link{aid}} and \code{\link{plot.ino}}.
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
#' @param a list contains inoculation grade data.
#' @return p-value signs indicating the significances between comparisons.
#' @seealso \code{\link{aid}}, and \code{\link{plot.ino}}.
#' @export
#' @examples
#' dat <- read.table(file, header = TRUE, check.names=FALSE)
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
    
    # assign p-value signs
    ino.p[i] <- t$p.value
    if (ino.p[i] < 0.005) {
      ino.ps[i] <- c("***")
    } else if (ino.p[i] < 0.01) {
      ino.ps[i] <- c("**")
    } else if (ino.p[i] < 0.05) {
      ino.ps[i] <- c("*")
    } else if (ino.p[i] < 0.1) {
      ino.ps[i] <- c("+")
    } else {
      ino.ps[i] <- c(" ")
    }
    
    i <- i + 1;

  }

  return(ino.ps)

}

#' Plot inoculation data
#'
#' Create a barplot for grade inoculation data
#' 
#' @param  a list contains inoculation grade data.
#' @seealso \code{\link{aid}} and \code{\link{ino.colors}}.
#' @export
#' @examples
#' dat <- read.table(file, header = TRUE, check.names=FALSE)
#' plot.ino(dat)
plot.ino <- function(ino) {

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
  legend("top", legend = colnames(ino), horiz = TRUE, fill=ino.colors(ng))

}

#' Analysis of inoculation data
#'
#' Analysis and illustration of inoculation data.
#' 
#' @param a text file contain inoculation data.
#' @return statistical analysis and illustration for inoculation data
#' @seealso \code{\link{grade.test}} and \code{\link{plot.ino}}.
#' @export
#' @examples
#' library(aid)
#' aid("~/ino.data")
aid <- function (file) {
  welcome()
  ino <- read.table(file, header = TRUE, check.names=FALSE)
  plot.ino(ino)
}
