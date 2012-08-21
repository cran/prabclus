\name{comp.test}
\alias{comp.test}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Compare species clustering and species groups}
\description{
  Tests for independence between a clustering and another grouping of species.
  This is simply an interface to \code{chisq.test}.
}

\usage{
comp.test(cl,spg)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{cl}{a vector of integers. Clustering of species (may be taken
    from \code{prabclust}).}
  \item{spg}{a vector of integers of the same length, groups of species.}
}
\details{
  \code{chisq.test} with simulated p-value is used.
}
\value{
  Output of \code{chisq.test}.
}

\references{
Hausdorf, B. and Hennig, C. (2003)  Biotic Element Analysis in
Biogeography. \emph{Systematic Biology} 52, 717-723.
}
\author{Christian Hennig
  \email{chrish@stats.ucl.ac.uk}
  \url{http://www.homepages.ucl.ac.uk/~ucakche}}


\seealso{
  \code{\link{chisq.test}}, \code{\link{prabclust}}.
}


\examples{
set.seed(1234)
g1 <- c(rep(1,34),rep(2,12),rep(3,15))
g2 <- sample(3,61,replace=TRUE)
comp.test(g1,g2)
}
\keyword{htest}% at least one, from doc/KEYWORDS

