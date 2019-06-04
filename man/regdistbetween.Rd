\name{regdistbetween}
\alias{regdistbetween}
\alias{print.regdistbetween}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Testing equality of within-groups and between-groups distances regression}
\description{
  Jackknife-based test for equality of two regressions between
  distances. Given two groups of objects, this tests whether the
  regression involving all distances is compatible with the regression
  involving within-group distances only.
}

\usage{
regdistbetween(dmx,dmy,grouping,groups=levels(as.factor(grouping))[1:2])

\method{print}{regdistbetween}(x,...)
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
  \item{x}{object of class \code{"regdistbetween"}.}
  \item{...}{optional arguments for print method.}
}

\details{
  The null hypothesis that the regressions based on all distances and
  based on within-group distances only are
  equal is tested using jackknife pseudovalues. This assumes that a
  single regression is appropriate at least for the within-group
  distances alone. The test statistic is the difference between fitted
  values with x (explanatory variable) fixed at the center of the
  between-group distances. The test is run one-sided, i.e., the null
  hypothesis is only rejected if the between-group distances are larger
  than expected under the null hypothesis, see below.   

  The test cannot be run in case that
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
  within-group distances. This is only rejected if the between-group
  distances are larger than expected under equality of regressions,
  because if they are smaller, this is not an indication against the
  groups belonging together genetically.

  If a joint regression on
  within-group distances is rejected by \code{regeqdist},
  \code{regdistbetweenone} can be
  used to test whether the between-group distances are at least
  compatible with the within-group distances of one of the groups, which
  can still be the case within a single species, see Hausdorf and Hennig (2019).
}


\value{
  list of class \code{"regdistbetween"} with components
  \item{pval}{p-value.} 
  \item{coeffdiff}{difference between regression fits (all distances
    minus within-group distances only) at \code{xcenterbetween}, see below.}  
  \item{condition}{condition numbers of regressions, see \code{\link{kappa}}.}
  \item{lmfit}{list. Output objects of \code{\link{lm}} within the two groups.}
  \item{jr}{output object of \code{\link[bootstrap]{jackknife}} for difference
    between regression fitted values at \code{xcenterbetween}.}
  \item{xcenter}{mean of within-groups distances of explanatory
    variable, used for centering.}
  \item{xcenterbetween}{mean of between-groups distances of explanatory
    variable (after centering by \code{xcenter}); at this point
    regression fitted values are computed.}
  \item{tstat}{t-statistic.}
  \item{tdf}{degrees of freedom of t-statistic.}
  \item{jackest}{jackknife-estimator of difference between regression
    fitted values at \code{xcenterbetween}.} 
  \item{jackse}{jackknife-standard error for
    \code{jackest}.}
  \item{jackpseudo}{vector of jacknife pseudovalues on which the test is
    based.}
  \item{testname}{title to be printed out when using
    \code{print.regdistbetween}.}
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
   \code{\link{regeqdist}}, \code{\link{regdistbetweenone}}
}

\examples{
  options(digits=4)
  data(veronica)
  ver.geo <- coord2dist(coordmatrix=veronica.coord[173:207,],file.format="decimal2")
  vei <- prabinit(prabmatrix=veronica[173:207,],distance="jaccard")
  loggeo <- log(ver.geo+quantile(as.vector(as.dist(ver.geo)),0.25))

  species <-c(rep(1,13),rep(2,22))
 
  rtest2 <-
  regdistbetween(dmx=loggeo,dmy=vei$distmat,grouping=species,groups=c(1,2))
  print(rtest2)
}
\keyword{htest}% __ONLY ONE__ keyword per line
\keyword{regression}% __ONLY ONE__ keyword per line
\keyword{spatial}% __ONLY ONE__ keyword per line
