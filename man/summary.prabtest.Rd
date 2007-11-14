\name{summary.prabtest}
\alias{summary.prabtest}
\alias{print.summary.prabtest}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Print and summary method for prabtest}
\description{
  Useful compression of the results of \code{prabtest}.
}
\usage{
\method{summary}{prabtest}(object, above.p=object$teststat \%in\%
         c("groups","inclusions","mean"),
         group.outmean=FALSE,...)

\method{print}{summary.prabtest}(x, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{object}{object of class \code{prabtest}.}
  \item{above.p}{logical. \code{TRUE} means that for output from
    \code{abundtest} the p-value is
    \code{p.above}, otherwise \code{p.below}.}
  \item{group.outmean}{logical. If \code{TRUE} and
    \code{object$teststat="groups"}, statistics concerning the mean of
    all dissimilarities are given out by \code{print.summary.prabtest}.}    
  \item{x}{object of class \code{summary.prabtest}.}
  \item{\dots}{no meaning, necessary for print and summary methods.}
}

\value{
  \code{summary.prabtest} produces an object of class
  \code{summary.prabtest}, which is a list with components
  \item{rrange}{range of the simulation results (test statistic values)
    of \code{object}.}
  \item{rmean}{mean of the simulation results (test statistic values)
    of \code{object}.}
  \item{datac, p.value, pd, tuning, teststat, distance, times, pdfnb,
    abund, sarlambda}{directly
    taken from \code{object}, see \code{prabtest} and \code{abundtest}.}
  \item{groupinfo}{if \code{object$teststat="groups"}, components
    \code{rrangeg} (matrix of group-wise ranges of test statistic
    value), \code{rmeang} (vector of group-wise means of test statistic
    value), \code{rrangem} (range over simulations of overall mean of
    within-group dissimilarities), \code{rmeanm} (mean over simulations
    of overall mean of within-group dissimilarities) are added to the
    list \code{object$groupinfo}, and this is given out.}
}

\author{Christian Hennig
  \email{chrish@stats.ucl.ac.uk}
  \url{http://www.homepages.ucl.ac.uk/~ucakche}}

\seealso{\code{\link{prabtest}}, \code{\link{abundtest}}}

\examples{
data(kykladspecreg)
data(nb)
set.seed(5678)
kprab <- prabinit(prabmatrix=kykladspecreg, neighborhood=nb)
kpt <- prabtest(kprab, times=5, pd=0.35, teststat="lcomponent")
# The parameters times and
# pd are only given to fasten the example execution.
summary(kpt)
}
\keyword{cluster}% at least one, from doc/KEYWORDS
\keyword{spatial}% __ONLY ONE__ keyword per line


