\name{initial.RndEM.worker}
\alias{initial.RndEM.worker}
\title{Initialization of RndEM for X.worker}
\description{
  This function implements RndEM procedure for model-based clustering
  based on \code{\link{X.worker}}.
}
\usage{
  initial.RndEM.worker(PARAM)
}
\arguments{
  \item{PARAM}{an original set of parameters generated
               by \code{\link{set.global}}.}
}
\details{
  The RndEM procedure is implemented by randomly picking
  \code{CONTROL$RndEM.iter} starting points from data \code{\link{X.worker}}
  and run one E-step to obtain the log likelihood.
  Then pick the starting point with the highest log likelihood as the
  best choice to pursue the MLEs in further EM iterations.

  This function repeatedly run \code{\link{initial.em.worker}} by
  \code{CONTROL$RndEM.iter} random starts and pick the best initializations
  from the random starts.
}
\value{
  The best initial starting points will be returned among all random
  starting points. The number of random starting points is assigned by
  \code{\link{set.global}} to a list variable \code{\link{CONTROL}}.
  See the help page
  of \code{\link{initial.em.worker}} and \code{\link{set.global}}
  for details.
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
  \code{\link{initial.em.worker}},
  \code{\link{em.step.worker}},
  \code{\link{aecm.step.worker}}, \cr
  \code{\link{apecm1.step.worker}},
  \code{\link{apecm2.step.worker}}.
}
\examples{
\dontrun{
# Examples can be found in the help pages of em.step.worker(),
# aecm.step.worker(), apecm1.step.worker(), and apecm2.step.worker().
}
}
\keyword{initialization}

