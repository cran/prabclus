\name{regeqdist}
\alias{regeqdist}
\alias{print.regeqdist}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Testing equality of two distance-regressions}
\description{
  Jackknife-based test for equality of two regressions between
  distance matrices.
}

\usage{
regeqdist(dmx,dmy,grouping,groups=levels(as.factor(grouping))[1:2])

\method{print}{regeqdist}(x,...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{dmx}{dissimilarity matrix or object of class
    \code{dist}. Explanatory dissimilarities (often these will be proper
    distances, but more general dissimilarities that do not
    necessarily fulfill the triangle inequality can be used, same for
    \code{dmy}).}
  \item{dmy}{dissimilarity matrix or object of class
    \code{dist}. Response dissimilarities.}
  \item{grouping}{something that can be coerced into a factor,
    defining the grouping of
    objects represented by the dissimilarities \code{dmx} and \code{dmy}
    (i.e., if \code{grouping} has length n, \code{dmx} and \code{dmy}
    must be dissimilarities between \code{n} objects).}
  \item{groups}{Vector of two, indicating the two groups defining the
    regressions to be compared in the test. These can be
    factor levels, integer numbers, or strings, depending on the entries
    of \code{grouping}.}
  \item{x}{object of class \code{"regeqdist"}.}
  \item{...}{optional arguments for print method.}
}

\details{
  The null hypothesis that the regressions within the two groups are
  equal is tested using jackknife pseudovalues independently in both
  groups allowing for potentially different variances of the
  pseudovalues, and aggregating as in Welch's t-test. Tests are run
  separately for intercept and slope and aggregated by Bonferroni's
  rule.

  The test cannot be run and many components will be \code{NA} in case that
  within-group regressions or jackknifed within-group regressions are
  ill-conditioned.

  This was implemented having in mind an application in which the
  explanatory distances represent geographical distances, the response
  distances are genetic distances, and groups represent species or
  species-candidates. In this application, for testing whether the
  regression patterns are compatble with the two groups behaving like a
  single species, one would first use \code{regeqdist} to test whether a
  joint regression for the within-group distances of both groups makes
  sense. If this is not rejected, \code{regdistbetween} is run to see
  whether the between-group distances are compatible with the
  within-group distances. On the other hand, if a joint regression on
  within-group distances is rejected, \code{regdistbetweenone} can be
  used to test whether the between-group distances are at least
  compatible with the within-group distances of one of the groups, which
  can still be the case within a single species, see Hausdorf and Hennig (2019).
}

\value{
  list of class \code{"regeqdist"} with components
  \item{pval}{p-values for intercept and slope.} 
  \item{coeffdiff}{vector of differences between groups (first minus
  second) for intercept and slope.}  
  \item{condition}{condition numbers of regressions, see \code{\link{kappa}}.}
  \item{lmfit}{list. Output objects of \code{\link{lm}} within the two groups.}
  \item{jr}{list of two lists of two; output object of
    \code{\link[bootstrap]{jackknife}} within the two groups for intercept and
    slope.}
  \item{xcenter}{mean of \code{dmx} within the two groups used for centering.}
  \item{tstat}{t-statistic.}
  \item{tdf}{vector of degrees of freedom of t-statistic according to
    Welch-Sattertwaithe approximation for intercept and slope.}
  \item{jackest}{jackknife-estimator of difference between regressions; 
    vector with intercept and slope difference.}
  \item{jackse}{vector with jackknife-standard errors for
    \code{jackest}, intercept and slope.}
  \item{jackpseudo}{list of two lists of vectors; jacknife pseudovalues
    within both groups for intercept and slope estimators.} 
  \item{groups}{see above.}
}
\references{
  Hausdorf, B. and Hennig, C. (2019) Species delimitation and
  geography. Submitted.

}

\author{Christian Hennig
  \email{christian.hennig@unibo.it}
  \url{https://www.unibo.it/sitoweb/christian.hennig/en}}

\seealso{
  \code{\link{regdistbetween}}, \code{\link{regdistbetweenone}}
}

\examples{
  options(digits=4)
  data(veronica)
  ver.geo <- coord2dist(coordmatrix=veronica.coord[173:207,],file.format="decimal2")
  vei <- prabinit(prabmatrix=veronica[173:207,],distance="jaccard")
  loggeo <- log(ver.geo+quantile(as.vector(as.dist(ver.geo)),0.25))

  species <-c(rep(1,13),rep(2,22))
  rtest <- regeqdist(dmx=loggeo,dmy=vei$distmat,grouping=species,groups=c(1,2))
  print(rtest)
}
\keyword{htest}% __ONLY ONE__ keyword per line
\keyword{regression}% __ONLY ONE__ keyword per line
\keyword{spatial}% __ONLY ONE__ keyword per line
