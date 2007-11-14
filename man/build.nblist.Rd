\name{build.nblist}
\alias{build.nblist}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Generate spatial weights from prabclus neighborhood list}
\description{
  This generates a \code{listw}-object as needed for estimation of a
  simultaneous autoregression model in package \code{spdep} from a
  neighborhood list of the type generated in \code{prabinit}.  
}
\usage{
build.nblist(prabobj,prab01=NULL,style="C")
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{prabobj}{object of class \code{prab}.}
  \item{prab01}{presence-absence matrix of same dimensions than the
  abundance matrix of \code{prabobj}. This specifies the presences and
  absences on which the presence/absence step of abundance-based tests
  is based (see details). If \code{NULL} (which is usually the only
  reasonable choice), \code{prab01} is computed in order to indicate
  the nonzeroes of \code{prabobj$prab}.}
  \item{style}{can take values "W", "B", "C", "U", and "S" though tests
    suggest that "C" should be chosen. See \code{\link[spdep]{nb2listw}}.} 
}
\value{
  A 'listw' object with the following members: 
\item{style}{see above.}
\item{neighbours}{the neighbours list in \code{spdep}-format.}
\item{weights}{the weights for the neighbours and chosen style, with
          attributes set to report the type of relationships (binary or
          general, if general the form of the glist argument), and
          style as above.} 
}

\author{Christian Hennig
  \email{chrish@stats.ucl.ac.uk}
  \url{http://www.homepages.ucl.ac.uk/~ucakche}}

\seealso{
  \code{\link[spdep]{nb2listw}} (which is called)
}

\examples{
# Not run; requires package spdep
# data(siskiyou)
# x <- prabinit(prabmatrix=siskiyou, neighborhood=siskiyou.nb,
#             distance="logkulczynski")
# build.nblist(x)
}
\keyword{spatial}% at least one, from doc/KEYWORDS

