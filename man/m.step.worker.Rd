\name{m.step.worker}
\alias{m.step.worker}
\title{Compute one M-step at current posterior probabilities}
\description{
  This function will perform one M-step based on current posterior
  probabilities.
  This is a core function of \code{\link{em.onestep.worker}}.
}
\usage{
  m.step.worker(PARAM)
}
\arguments{
  \item{PARAM}{a set of parameters.}
}
\details{
  This function will base on the current posterior probabilities
  \code{\link{Z.worker}} to estimate the parameters \code{\link{PARAM}}
  mainly including mixing proportions \code{\link{ETA}},
  centers of clusters \code{\link{MU}}, and
  dispersions of clusters \code{\link{SIGMA}}.
}
\value{
  Returning a new \code{\link{PARAM}} which maximizes the complete data log
  likelihood for the current iteration.
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
  \code{\link{em.onestep.worker}},
  \code{\link{e.step.worker}}.
}
\examples{
\dontrun{
# This is a core function for em.onestep.worker()
# see the source code for details.
}
}
\keyword{core function}

