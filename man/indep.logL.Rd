\name{indep.logL}
\alias{indep.logL}
\title{Independent function for log likelihood}
\description{
  This function is for debugging only and for checking if the observed data
  log likelihood is consistent for each EM iteration.
}
\usage{
  indep.logL(PARAM)
}
\arguments{
  \item{PARAM}{a set of parameters.}
}
\details{
  This function will provide an observed data log likelihood based on the
  current parameter \code{\link{PARAM}}.
  This function will take in information from global, but
  no global variables will be updated by this function.

  This function also don't take care the numerical issues, so the return
  value may be inaccurate sometimes.
}
\value{
  An observed data log likelihood will be returned.
  This value can quickly compare with the log likelihood computed inside
  \code{\link{em.onestep.worker}}. Small difference is allowed, but
  large difference indicates bugs of code or illness of data.
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
  \code{\link{em.onestep.worker}}.
}
\examples{
\dontrun{
# This is a core function for em.estep.worker()
# see the source code for details.
# Reset CONTROL$debug to turn on this function
# automatically for each EM iteration.
}
}
\keyword{debugging function}

