\name{initial.center.worker}
\alias{initial.center.worker}
\title{Initialization for K-means algorithm}
\description{
  This function implements simple initialization for K-means algorithm
  by randomly picking samples from data based on X.worker.
}
\usage{
  initial.center.worker(PARAM, MU = NULL)
}
\arguments{
  \item{PARAM}{an original set of parameters generated
               by \code{\link{set.global}}.}
  \item{MU}{a center matrix with dim = \eqn{p \times K}{\code{p * K}}.}
}
\details{
  If \code{\link{MU}} is given, then the center will be assigned according.
}
\value{
  An initial set of parameters \code{\link{PARAM}} will be returned. 
}
\references{
  High Performance Statistical Computing Website:
  \url{http://thirteen-01.stat.iastate.edu/snoweye/hpsc/}
}
\author{
  Wei-Chen Chen \email{wccsnow@gmail.com} and George Ostrouchov.
}
\seealso{
  \code{\link{set.global}},
  \code{\link{initial.RndEM.worker}}.
}
\examples{
\dontrun{
# Examples can be found in the help page of kmeans.step.worker().
}
}
\keyword{initialization}

