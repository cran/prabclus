\name{specgroups}
\alias{specgroups}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Average within-group distances for given groups}
\description{
  Generates average within-group distances (overall and group-wise) from
  a dissimilarity matrix and a given grouping.  
}
\usage{
  specgroups(distmat,groupvector, groupinfo) 
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{distmat}{dissimilarity matrix or \code{dist}-object.}
  \item{groupvector}{integer vector. For every row of \code{distmat},
  a number indicating the group membership.}
  \item{groupinfo}{list with components \code{lg} (levels of
    \code{groupvector}), \code{ng} (number of groups), \code{nsg}
    (vector of group sizes).}
}
\value{
  A list with parameters 
\item{overall}{overall average within-groups dissimilarity.}
\item{gr}{vector of group-wise average within-group dissimilarities.}
}

\author{Christian Hennig
  \email{chrish@stats.ucl.ac.uk}
  \url{http://www.homepages.ucl.ac.uk/~ucakche}}

\examples{
  data(siskiyou)
  x <- prabinit(prabmatrix=siskiyou, neighborhood=siskiyou.nb,
              distance="logkulczynski")  
  groupvector <- as.factor(siskiyou.groups)
  ng <- length(levels(groupvector))
  lg <- levels(groupvector)
  nsg <- numeric(0)
    for (i in 1:ng) nsg[i] <- sum(groupvector==lg[i])
  groupinfo <- list(lg=lg,ng=ng,nsg=nsg)
  specgroups(x$distmat,groupvector,groupinfo)
}
\keyword{cluster}% at least one, from doc/KEYWORDS

