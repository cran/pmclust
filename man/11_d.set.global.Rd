\name{.1.set.global}
\alias{set.global}
\alias{X.spmd}
\title{Set Global Variables According to X.spmd}
\description{
  This function will set several sets of variables globally according to
  the \code{X.spmd}.
}
\usage{
  set.global(K = 2, PARAM = NULL,
      method = c("em", "aecm", "apecm1", "apecm2", "kmeans"),
      RndEM.iter = 10)
}
\arguments{
  \item{K}{an original set of parameters generated
           by \code{set.global}.}
  \item{PARAM}{an original set of parameters generated
               by \code{set.global}.}
  \item{method}{an original set of parameters generated
                by \code{set.global}.}
  \item{RndEM.iter}{number of RndEM iterations.}
}
\details{
  WARNING: a global variable named \code{X.spmd} should be set before
  calling this function where \code{X.spmd} is a matrix containing data
  with dimension \code{N.spmd * p}. i.e. \code{N.spmd} observations and
  \code{p} variables.
}
\value{
  A new set of \code{\link{PARAM}} will be returned and several global
  variables will be set according to the data \code{X.spmd}.

  Use \code{\link{readme}} to see all global variables set by this function.
}
\references{
  High Performance Statistical Computing Website:
  \url{http://thirteen-01.stat.iastate.edu/snoweye/hpsc/}
}
\author{
  Wei-Chen Chen \email{wccsnow@gmail.com} and George Ostrouchov.
}
\seealso{
  \code{\link{em.step.spmd}},
  \code{\link{aecm.step.spmd}},
  \code{\link{apecm1.step.spmd}}, \cr
  \code{\link{apecm2.step.spmd}}.
}
\examples{
\dontrun{
# Examples can be found in the help pages of em.step.spmd(),
# aecm.step.spmd(), apecm1.step.spmd(), and apecm2.step.spmd().
}
}
\keyword{global variables}

