\name{initial.em.worker}
\alias{initial.em.worker}
\title{Initialization for model-based clustering}
\description{
  This function implements one E-step for one random initial start.
}
\usage{
  initial.em.worker(PARAM, MU = NULL)
}
\arguments{
  \item{PARAM}{an original set of parameters generated
               by \code{set.global}.}
  \item{MU}{a center matrix with dim = \eqn{p \times K}{\code{p * K}}.}
}
\details{
  This function takes \code{\link{X.worker}} from the global environment and
  randomly pick \eqn{K}{K} of them as the centers of \eqn{K}{K} groups.
  If \code{\link{MU}} is specified, then this \code{\link{MU}}
  will be the centers.
  The default identity dispersion in \code{PARAM$SIGMA} will be used.
  Then, one E-step will be called to obtain the log likelihood and new
  classification will be updated.

  This function is used to implement the RndEM procedure for more
  elaborate initialization scheme in \code{\link{initial.RndEM.worker}}.
  Potentially, several random starts should be tried before running EM
  algorithms. This can benefit in two aspects including:
  shorter convergent iterations and better classification results.
}
\value{
  An initial set of parameters \code{\link{PARAM}} will be returned.
}
\references{
  High Performance Statistical Computing Website:
  \url{http://thirteen-01.stat.iastate.edu/snoweye/hpsc/}

  Maitra, R. (2009)
  \dQuote{Initializing partition-optimization algorithms},
  \emph{IEEE/ACM Transactions on Computational Biology and Bioinformatics},
  \bold{6:1}, 114-157.
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
# Examples can be found in the help pages of em.step.worker(),
# aecm.step.worker(), apecm1.step.worker(), and apecm2.step.worker().
}
}
\keyword{initialization}

