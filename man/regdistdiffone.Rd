\name{regdistdiffone}
\alias{regdistdiffone}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Regression difference within reference group and between-group dissimilarities}
\description{
  Given two dissimilarity matrices \code{dmx} and \code{dmy}, an indicator
  vector \code{x} and a grouping, this computes the difference between
  standard least squares regression predictions at point
  \code{xcenterbetween}. The regressions are based on the dissimilarities
  in \code{dmx} vs. \code{dmy} for objects indicated in
  \code{x}. \code{grouping} indicates the two groups, and the difference
  is computed between regressions based on (a) the within-group
  distances of the reference group \code{rgroup} and (b) these together
  with the between-group distances. 
}

\usage{

regdistdiffone(x,dmx,dmy,grouping,xcenter=0,xcenterbetween=0,rgroup)

}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{vector of logicals of length of the number of objects on which
    dissimilarities \code{dmx} and \code{dmy} are based.} 
  \item{dmx}{dissimilarity matrix or object of class
    \code{\link{dist}}. Explanatory  dissimilarities.}
  \item{dmy}{dissimilarity matrix or object of class
    \code{\link{dist}}. Response dissimilarities.}
  \item{grouping}{vector of length of the number of objects on which
    dissimilarities \code{dmx} and \code{dmy} are based. Grouping
    vector. Regressions will be based on the first two values that
    appear in \code{unique(grouping[x])} (note that objects that are not
    assigned to one of these groups will be ignored); normally
    \code{grouping} should indicate only
    two groups on the objects with \code{x=TRUE}, and then these are used.}
  \item{xcenter}{numeric. Dissimilarities \code{dmx} are centered by
    this, i.e., this value is subtracted from the dissimilarities before
    regression.}
  \item{xcenterbetween}{numeric. This specifies the x- (dissimilarity)
    value at which predictions from the two regressions are
    compared. Note that this is interpreted as after centering by
    \code{xcenter}.}
  \item{rgroup}{one of the values of \code{grouping}, specifying the
    reference group.}

}

\value{
  Difference between
  standard least squares regression predictions for the two regressions at point
  \code{xcenterbetween}.  
}

\references{
  Hausdorf, B. and Hennig, C. (2019) Species delimitation and
  geography. Submitted.

}

\seealso{
  \code{\link{regdistbetweenone}}
}
  
\author{Christian Hennig
  \email{christian.hennig@unibo.it}
  \url{https://www.unibo.it/sitoweb/christian.hennig/en}}



\examples{
  options(digits=4)
  data(veronica)
  ver.geo <- coord2dist(coordmatrix=veronica.coord[173:207,],
    file.format="decimal2")
  vei <- prabinit(prabmatrix=veronica[173:207,],distance="jaccard")

  species <-c(rep(1,13),rep(2,22))
  regdistdiffone(rep(TRUE,35),ver.geo,vei$distmat,grouping=species,
    xcenter=0,xcenterbetween=100,rgroup=2)

}
\keyword{regression}% __ONLY ONE__ keyword per line
\keyword{spatial}% __ONLY ONE__ keyword per line
