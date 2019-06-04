\name{regdist}
\alias{regdist}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Regression between subsets of dissimilarity matrices}
\description{
  Given two dissimilarity matrices \code{dmx} and \code{dmy} and an indicator
  vector \code{x}, this computes a standard least squares regression
  between the dissimilarity between objects indicated in \code{x}.
}

\usage{

regdist(x,dmx,dmy,xcenter=0,param)

}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{vector of logicals of length of the number of objects on which
    dissimilarities \code{dmx} and \code{dmy} are based.} 
  \item{dmx}{dissimilarity matrix or object of class
    \code{\link{dist}}. Explanatory  dissimilarities.}
  \item{dmy}{dissimilarity matrix or object of class
    \code{\link{dist}}. Response dissimilarities.}
  \item{xcenter}{numeric. Dissimilarities \code{dmx} are centered by
    this, i.e., this value is subtracted from the dissimilarities before
    regression.}
  \item{param}{1 or 2 or \code{NULL}. If 1 or 2, only the first or
    second parameter (intercept or slope) of the regression is given
    out.}
}

\value{
  If \code{param=NULL}, the output object of \code{\link{lm}}. If
  \code{param=1} the intercept. If
  \code{param=2} the slope. 
}
\references{
  Hausdorf, B. and Hennig, C. (2019) Species delimitation and
  geography. Submitted.

}

\author{Christian Hennig
  \email{christian.hennig@unibo.it}
  \url{https://www.unibo.it/sitoweb/christian.hennig/en}}



\examples{
  options(digits=4)
  data(veronica)
  ver.geo <- coord2dist(coordmatrix=veronica.coord[1:20,],file.format="decimal2")
  vei <- prabinit(prabmatrix=veronica[1:20,],distance="jaccard")
  regdist(c(rep(TRUE,10),rep(FALSE,10)),ver.geo,vei$distmat,param=1)
}
\keyword{regression}% __ONLY ONE__ keyword per line
\keyword{spatial}% __ONLY ONE__ keyword per line
