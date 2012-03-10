\name{kmeans.update.class.spmd}
\alias{kmeans.update.class.spmd}
\title{Update CLASS.spmd Based on the K-Means Final Iteration}
\description{
  Update \code{\link{CLASS.spmd}} based on the final iteration of K-means
  algorithms.
}
\usage{
  kmeans.update.class.spmd()
}
\details{
  This function takes \code{\link{Z.spmd}} from the global environment and
  update \code{\link{CLASS.spmd}}, and provides the identification of groups
  for all data.
}
\value{
  \code{\link{CLASS.spmd}} will be updated.
}
\references{
  High Performance Statistical Computing Website:
  \url{http://thirteen-01.stat.iastate.edu/snoweye/hpsc/}
}
\author{
  Wei-Chen Chen \email{wccsnow@gmail.com} and George Ostrouchov.
}
\seealso{
  \code{\link{em.update.class.spmd}}.
}
\examples{
\dontrun{
# Examples can be found in the help pages of kmeans.step.spmd().
}
}
\keyword{programming}

