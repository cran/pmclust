\name{e.step.worker}
\alias{e.step.worker}
\title{Compute one E-step and log likelihood at current parameters}
\description{
  This function will perform one E-step based on current parameters.
  This is a core function of \code{\link{em.onestep.worker}}.
}
\usage{
  e.step.worker(PARAM, update.logL = TRUE)
}
\arguments{
  \item{PARAM}{a set of parameters.}
  \item{update.logL}{TRUE for update observed data log likelihood.}
}
\details{
  This function will base on the current parameter to 
  compute the densities for all observations for all
  \code{K} components, and update the \code{\link{Z.worker}} matrix.
  If the \code{update.logL} is true, then the log likelihood
  \code{\link{W.worker.rowSums}} will be also updated before the end
  of this function.

  Sum of \code{\link{W.worker.rowSums}} of all processors will be the
  observed data log likelihood for the current iteration.
}
\value{
  Several global variables will be overwrote after this call including
  \code{\link{Z.worker}}, \code{\link{W.worker.rowSums}},
  \code{\link{W.worker}}, \code{\link{U.worker}},
  and \code{\link{Z.colSums}}.
}
\section{Computing Issues}{
  Since the clusters can be degenerated or highly flat, these cause very
  large positive or negative exponents in densities.
  The log likelihood will tend to be inaccurate (not finite).
  Since the mixture structures can be over fit, this also cause very
  tiny mixing proportions.
  The poster probabilities can also unstable (NaN).

  These can be solved by rescaling the range of exponents carefully
  and adjust the scaling factor on the log values.
  See \code{\link{CONTROL}} for details about constrains on E- and M-steps.
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
  \code{\link{m.step.worker}}.
}
\examples{
\dontrun{
# This is a core function for em.onestep.worker()
# see the source code for details.
}
}
\keyword{core function}

