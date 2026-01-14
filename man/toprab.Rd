\name{toprab}
\alias{toprab}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Convert abundance matrix into presence/absence matrix}
\description{
  Converts abundance matrix into binary (logical) presence/absence
  matrix (\code{TRUE} if
  abundance>0). 
}
\usage{
toprab(prabobj)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{prabobj}{object of class \code{prab}.}
}
\value{
  Logical matrix with same dimensions as \code{prabobj$prab} as described above.
}

\author{Christian Hennig
  \email{christian.hennig@unibo.it}
  \url{https://www.unibo.it/sitoweb/christian.hennig/en}}


\examples{
  data(siskiyou)
  x <- prabinit(prabmatrix=siskiyou, neighborhood=siskiyou.nb,
             distance="none")
  toprab(x)
}
\keyword{manip}% at least one, from doc/KEYWORDS

