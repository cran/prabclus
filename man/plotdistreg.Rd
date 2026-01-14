\name{plotdistreg}
\alias{plotdistreg}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Plots for within-groups and between-groups distance regression}
\description{
  Visualisation of various regressions on distance (or dissimilarity)
  data where objects are from two groups.
}

\usage{
plotdistreg(dmx,dmy,grouping,groups=levels(as.factor(grouping))[1:2],
                        cols=c(1,2,3,4),
                        pchs=rep(1,3),
                        ltys=c(1,2,1,2),
                        individual=TRUE,jointwithin=TRUE,jointall=TRUE,
                        oneplusjoint=TRUE,jittering=TRUE,bcenterline=TRUE,
                        xlim=NULL,ylim=NULL,xlab="geographical distance",
                        ylab="genetic distance",...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dmx}{dissimilarity matrix or object of class
    \code{dist}. Explanatory dissimilarities (often these will be proper
    distances, but more general dissimilarities that do not
    necessarily fulfill the triangle inequality can be used, same for \code{dmy}).}
  \item{dmy}{dissimilarity matrix or object of class
    \code{dist}. Response dissimilarities.}
  \item{grouping}{something that can be coerced into a factor,
    defining the grouping of
    objects represented by the dissimilarities \code{dmx} and \code{dmy}
    (i.e., if \code{grouping} has length n, \code{dmx} and \code{dmy}
    must be dissimilarities between \code{n} objects).}
  \item{groups}{Vector of two levels. The two groups defining the
    regressions to be compared in the test. These can be
    factor levels, integer numbers, or strings, depending on the entries
    of \code{grouping}.}
  \item{cols}{vector of four colors (or color numbers) to be used for
    plotting distances 
    and regression lines within the first group, within the second group,
    distances between groups, and a line marking the center of the
    between-groups explanatory distances, see \code{col}-argument of
    \code{\link{par}}.}
  \item{pchs}{vector of three plot symbols (or numbers) to be used for
    plotting distances within the first group, within the second group,
    and distances between groups, see \code{pch}-argument of
    \code{\link{par}}.}
  \item{ltys}{vector of line type numbers to be used for single group
    within-group regression, both groups combined within-group
    regression, regression with all distances, and regression combining
    within-groups distances of one group with between-groups distances,
    see \code{lty}-argument of
    \code{\link{par}}.}
  \item{individual}{if \code{TRUE}, within-groups distances regression
    lines are shown for both groups.}
  \item{jointwithin}{if \code{TRUE}, the within-groups distances regression
    line for both groups combined is shown.}
  \item{jointall}{if \code{TRUE}, the regression
    line based on all distances is shown.}
  \item{oneplusjoint}{if \code{TRUE}, the regression lines combining
    within-groups distances of one group with between-groups distances
    are shown (the colors of these are the colors of the individual
    groups, the first two components of the \code{cols}-argument).}
  \item{jittering}{if \code{TRUE}, points are jittered to avoid
    overplotting.}
  \item{bcenterline}{if \code{TRUE}, a line is plotted to mark the
    center of the between-groups distances on the explanatory variable.}
  \item{xlim}{to be passed on to \code{\link{plot}}; the default is
    determined from the involved distances.}
  \item{ylim}{to be passed on to \code{\link{plot}}; the default is
    determined from the involved distances.}
  \item{xlab}{to be passed on to \code{\link{plot}}.}
  \item{ylab}{to be passed on to \code{\link{plot}}.}  
  \item{...}{optional arguments to be passed on to \code{\link{plot}}.}
}

\references{
  Hausdorf, B. and Hennig, C. (2019) Species delimitation and
  geography. Submitted.

}

\author{Christian Hennig
  \email{christian.hennig@unibo.it}
  \url{https://www.unibo.it/sitoweb/christian.hennig/en}}

\seealso{
  \code{\link{regeqdist}}, \code{\link{regdistbetween}},
  \code{\link{regdistbetweenone}}, \code{\link{regdistdiffone}}
}

\examples{
  options(digits=4)
  data(veronica)
  ver.geo <- coord2dist(coordmatrix=veronica.coord[173:207,],file.format="decimal2")
  vei <- prabinit(prabmatrix=veronica[173:207,],distance="jaccard")

  species <-c(rep(1,13),rep(2,22))
  loggeo <- log(ver.geo+quantile(as.vector(as.dist(ver.geo)),0.25))
  plotdistreg(dmx=loggeo,dmy=vei$distmat,grouping=species,
  jointwithin=FALSE,jointall=FALSE,groups=c(1,2))
  legend(5,0.75,c("within species 1",
  "within species 2","species 1 and between","species 2 and between"),lty=c(1,1,2,2),col=c(1,2,1,2))
  plotdistreg(dmx=loggeo,dmy=vei$distmat,grouping=species,
  jointwithin=TRUE,jointall=TRUE,oneplusjoint=FALSE,groups=c(1,2))
  legend(5,0.75,c("within species 1",
  "within species 2","all distances","all within species"),lty=c(1,1,1,2),col=c(1,2,3,3))
  

}
\keyword{htest}% __ONLY ONE__ keyword per line
\keyword{regression}% __ONLY ONE__ keyword per line
\keyword{spatial}% __ONLY ONE__ keyword per line
