\name{em.update.class.spmd}
\alias{em.update.class.spmd}
\title{Update CLASS.spmd Based on the Final Iteration}
\description{
  Update \code{\link{CLASS.spmd}} based on the final iteration of EM-like
  algorithms.
}
\usage{
  em.update.class.spmd()
}
\details{
  This function takes \code{\link{Z.spmd}} from the global environment
  \code{\link[base]{.GlobalEnv}} and update \code{\link{CLASS.spmd}},
  and provides the identification of groups for all data.
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
  \code{\link{kmeans.update.class.spmd}}.
}
\examples{
\dontrun{
# Examples can be found in the help pages of em.step.spmd(),
# aecm.step.spmd(), apecm1.step.spmd(), and apecm2.step.spmd().
}
}
\keyword{programming}

