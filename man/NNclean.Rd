\name{NNclean}
\alias{NNclean}
\alias{print.nnclean}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Nearest neighbor based clutter/noise detection}
\description{
  Detects if data points are noise or part of a cluster,
  based on a Poisson process model.
}
\usage{
NNclean(data, k, distances = NULL, edge.correct = FALSE, wrap = 0.1,
convergence = 0.001, plot=FALSE)

\method{print}{nnclean}(x, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{data}{numerical matrix or data frame.}
  \item{k}{integer. Nmber of considered nearest neighbors per point.}
  \item{distances}{distance matrix object of class \code{dist}. If
    specified, it is used instead of computing distances from the data.}
  \item{edge.correct}{logical. If \code{TRUE} and the data is
    two-dimensional, neighbors for points at the edges of the parent
    region of the noise Poisson process are determined after wrapping
    the region onto a toroid.}
  \item{wrap}{numerical. If \code{edge.correct=TRUE}, points in a
    strip of size \code{wrap*range} along the edge for each variable
    are candidates for
    being neighbors of points from the opposite.}
  \item{convergence}{numerical. Convergence criterion for EM-algorithm.}
  \item{plot}{logical. If \code{TRUE}, a histogram of the distance to
    kth nearest neighbor and fit is plotted.}
  \item{x}{object of class \code{nnclean}.}
  \item{...}{necessary for print methods.}
}
\details{
  The assumption is that the noise is distributed as a homogeneous
  Poisson process  on a certain region and the clusters are distributed
  as a homogeneous Poisson process with larger intensity on a
  subregion (disconnected in case of more than one cluster).
  The distances are then distributed according to a mixture of two
  transformed Gamma distributions, and this mixture is estimated via the
  EM-algorithm. The points are assigned to noise or cluster component
  by use of the estimated a posteriori probabilities. 
}
\value{
  \code{NNclean} returns a list of class \code{nnclean} with components
  \item{z}{0-1-vector of lenght of the number of data points. 1 means
    cluster, 0 means noise.}
  \item{probs}{vector of estimated a priori probabilities for each point
    to belong to the cluster component.}
  \item{k}{see above.}
  \item{lambda1}{intensity parameter of cluster component.}
  \item{lambda2}{intensity parameter of noise component.}
  \item{p}{estimated probability of cluster component.}
  \item{kthNND}{distance to kth nearest neighbor.}
}
\references{
Byers, S. and Raftery, A. E. (1998) Nearest-Neighbor Clutter
Removal for Estimating Features in Spatial Point Processes,
\emph{Journal of the American Statistical Association}, 93, 577-584.
}
\author{R-port by Christian Hennig
  \email{hennig@math.uni-hamburg.de}
  \url{http://www.math.uni-hamburg.de/home/hennig/},\cr
  original Splus package by S. Byers and A. E. Raftery.
}

\note{The software can be freely used for non-commercial purposes, and can
be freely distributed for non-commercial purposes only.} 

\examples{
library(mclust)
data(chevron)
nnc <-  NNclean(chevron[,2:3],15,plot=TRUE)
plot(chevron[,2:3],col=1+nnc$z)
}

\keyword{multivariate}% at least one, from doc/KEYWORDS
\keyword{cluster}% __ONLY ONE__ keyword per line
