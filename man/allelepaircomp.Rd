\name{allelepaircomp}
\alias{allelepaircomp}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Internal: compares two pairs of alleles}
\description{
  Used for computation of the genetic distances \code{\link{alleledist}}.
%-  \code{\link{neidist}}, \code{\link{chorddist}}.
}
\usage{
  allelepaircomp(allelepair1,allelepair2,method="sum")
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{allelepair1}{vector of two allele codes (usually characters), or
    \code{NA}.}
  \item{allelepair2}{vector of two allele codes (usually characters), or
    \code{NA}.}
  \item{method}{one of \code{"sum"} or \code{"geometrical"}.}
}

\value{
  If \code{method=="sum"}, number of shared alleles (0, 1 or 2), or
  \code{NA}. If \code{method=="geometrical"}, 0, 0.5, \code{sqrt(0.5)}
  (in case that one of the allelepairs is double such as in
  \code{c("A","B"),c("A","A")}) or 1, or
  \code{NA}. 
}

\author{Christian Hennig
  \email{chrish@stats.ucl.ac.uk}
  \url{http://www.homepages.ucl.ac.uk/~ucakche}}

\seealso{
  \code{\link{alleledist}}
%-  , \code{\link{neidist}}, \code{\link{chorddist}}.
}

\examples{
  allelepaircomp(c("A","B"),c("A","C"))
}
\keyword{cluster}% at least one, from doc/KEYWORDS
