\name{nastats}
\alias{nastats}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Missing values statistics for matrix}
\description{
  Computes column-wise and row-wise numbers of missing values. 
}

\usage{
  nastats(amatrix, nastr="--")
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{amatrix}{(any) matrix.}
  \item{nastr}{missing value indicator.}
}
  
\value{
  A list with components
  \item{narow}{vector of row-wise numbers of mixxing values.}
  \item{nacol}{vector of column-wise numbers of mixxing values.}
}

\author{Christian Hennig
  \email{christian.hennig@unibo.it}
  \url{https://www.unibo.it/sitoweb/christian.hennig/en}}

\examples{
  xx <- cbind(c(1,2,3),c(0,0,1),c(5,3,1))
  nastats(xx,nastr=0)
}
\keyword{manip}% at least one, from doc/KEYWORDS




