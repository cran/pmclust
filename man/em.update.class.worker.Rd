\name{em.update.class.worker}
\alias{em.update.class.worker}
\title{Update CLASS.worker based on the final iteration}
\description{
  Update \code{\link{CLASS.worker}} based on the final iteration of EM-like
  algorithms.
}
\usage{
  em.update.class.worker()
}
\details{
  This function takes \code{\link{Z.worker}} from the global environment
  \code{\link[base]{.GlobalEnv}} and update \code{\link{CLASS.worker}},
  and provides the identification of groups for all data.
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
  \code{\link{kmeans.update.class.worker}}.
}
\examples{
\dontrun{
# Examples can be found in the help pages of em.step.worker(),
# aecm.step.worker(), apecm1.step.worker(), and apecm2.step.worker().
}
}
\keyword{programming}

