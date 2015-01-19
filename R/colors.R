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