\name{build.charmatrix}
\alias{build.charmatrix}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Internal: create character matrix out of allele list}
\description{
  For use in \code{\link{alleleinit}}.
  Creates a matrix of characters in which there are two rows for
    every individual corresponding to the two alleles in every locus
    (column) out of a list of lists, such as required by
    \code{\link{alleledist}}.
}

\usage{
   build.charmatrix(allelelist,n.individuals,n.variables)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{allelelist}{A list of lists. In the "outer" list, there are
  \code{n.variables} lists, one for each locus. In the "inner" list, for every
    individual there is a vector of two codes (typically characters, see
    \code{\link{alleleinit}}) for the two alleles in that locus.}  
  \item{n.individuals}{integer. Number of individuals.}
  \item{n.variables}{integer. Number of loci.}
}
  
\value{
  A matrix of characters in which there are two rows for
    every individual corresponding to the two alleles in every locus
    (column). 
}

\author{Christian Hennig
  \email{chrish@stats.ucl.ac.uk}
  \url{http://www.homepages.ucl.ac.uk/~ucakche}}

\seealso{\code{\link{alleleinit}}, \code{\link{unbuild.charmatrix}}}

\examples{
  alist <- list()
  alist[[1]] <- list(c("A","A"),c("B","A"),c(NA,NA))
  alist[[2]] <- list(c("A","C"),c("B","B"),c("A","D"))
  build.charmatrix(alist,3,2)
}
\keyword{manip}% at least one, from doc/KEYWORDS




