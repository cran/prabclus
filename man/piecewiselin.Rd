\name{piecewiselin}
\alias{piecewiselin}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Piecewise linear transformation for distance matrices}
\description{
  Piecewise linear transformation for distance matrices, utility
  function for \code{geco}.
}
\usage{
piecewiselin(distmatrix, maxdist=0.1*max(distmatrix))
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{distmatrix}{symmetric (non-negative) distance matrix.}
  \item{maxdist}{non-negative numeric. Larger distances are transformed
    to constant 1.}
}
\details{
  Transforms large distances to 1, 0 to 0 and continuously linear between 0 and
  \code{maxdist}.
}
\value{
  A symmetrical matrix.
}
  
\author{Christian Hennig
  \email{chrish@stats.ucl.ac.uk}
  \url{http://www.homepages.ucl.ac.uk/~ucakche}}


\seealso{
  \code{\link{geco}}
}

\examples{
data(waterdist)
piecewiselin(waterdist)
}
\keyword{cluster}% at least one, from doc/KEYWORDS
\keyword{spatial}% __ONLY ONE__ keyword per line

