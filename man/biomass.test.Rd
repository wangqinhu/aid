% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/test.R
\name{biomass.test}
\alias{biomass.test}
\title{Statistical test for biomass data}
\usage{
biomass.test(bio, alternative, paired)
}
\arguments{
\item{alternative}{a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less". You can specify just the initial letter.}

\item{paired}{a character string specifying if the experiment is paired, can be TRUE or FALSE (default). You can specify just the initial letter.}

\item{ino}{a list contains inoculation biomass data.}
}
\value{
p-value symbols indicating the significances between comparisons.
}
\description{
Perform t test on biomass (qPCR) data.
}
\examples{
bio<-c(0.82, 3.14, 0.88, 3.21, 0.85, 3.20)
dim(bio)<-c(2,3)
biomass.test(bio, alternative = "t", paired = FALSE)
}
\seealso{
\code{\link{aid}}, \code{\link{biomass.barplot}}, \code{\link{lesion.test}} and \code{\link{grade.test}}.
}

