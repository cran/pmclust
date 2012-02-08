\name{kmeans.update.class.worker}
\alias{kmeans.update.class.worker}
\title{Update CLASS.worker based on the K-means final iteration}
\description{
  Update \code{\link{CLASS.worker}} based on the final iteration of K-means
  algorithms.
}
\usage{
  kmeans.update.class.worker()
}
\details{
  This function takes \code{\link{Z.worker}} from the global environment and
  update \code{\link{CLASS.worker}}, and provides the identification of groups
  for all data.
}
\value{
  \code{\link{CLASS.worker}} will be updated.
}
\references{
  High Performance Statistical Computing Website:
  \url{http://thirteen-01.stat.iastate.edu/snoweye/hpsc/}
}
\author{
  Wei-Chen Chen \email{wccsnow@gmail.com} and George Ostrouchov.
}
\seealso{
  \code{\link{em.update.class.worker}}.
}
\examples{
\dontrun{
# Examples can be found in the help pages of kmeans.step.worker().
}
}
\keyword{programming}

