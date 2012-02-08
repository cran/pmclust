\name{generate.basic.worker}
\alias{generate.basic.worker}
\title{Generate basic examples for testing}
\description{
  This function will generate a small set of data for testing algorithms.
}
\usage{
  generate.basic.worker(N.allworkers, N.worker, N.K.worker,
                        N, p, K, seed)
}
\arguments{
  \item{N.allworkers}{a collection of sample sizes for all
                      \eqn{S}{S} processors, i.e. a vector of
                      length \eqn{S}{S}.}
  \item{N.worker}{total sample size of given processor.}
  \item{N.K.worker}{sample size of each clusters given processor, i.e.
                    sum over \code{N.K.worker} is \code{N.worker},
                    a vector of length
                    \eqn{K}{K}.}
  \item{N}{total sample size across all \eqn{S}{S} processors, i.e.
           sum over \code{N.worker} is \code{N}.}
  \item{p}{dimension of data \code{\link{X.worker}},
           i.e. \code{ncol(\link{X.worker})}.}
  \item{K}{number of clusters.}
  \item{seed}{a seed for random numbers and processor's rank will be added
              later.}
}
\details{
  For all \eqn{S}{S} processors, this function will generate in total
  \eqn{N}{N} observations from \eqn{K}{K} clusters in \eqn{p}{p} dimensions.

  The clusters centers and dispersions are generated automatically inside
  the code. Currently, it is not allowed for users to change, but it is
  not difficult to specify them by mimicking this code.
}
\value{
  A set of simulated data and information will be returned in a list variable
  including:
\tabular{ll}{\code{K}               \tab number of clusters, as the input \cr
             \code{p}               \tab dimension of data
                                         \code{\link{X.worker}},
                                         as the input \cr
             \code{N}               \tab total sample size, as the input \cr
             \code{N.allworkers}    \tab a collection of sample sizes for all
                                         \eqn{S}{S} processors, as the input \cr
             \code{N.worker}        \tab total sample size of given processor,
                                         as the input \cr
             \code{N.K.worker}      \tab sample size of each clusters given
                                         processor, as the input \cr
             \code{seed}            \tab a seed for random numbers, as the
                                         input \cr
             \code{\link{X.worker}} \tab generated data set with dimension with
                                         dimension \code{N.worker * p} \cr
             \code{\link{CLASS.worker}}
                                    \tab true id of each data, a vector of
                                         length \code{N.worker}
                                         and has values from 1 to \code{K} \cr
             \code{N.CLASS.worker}  \tab true sample size of each clusters, a
                                         vector of length \code{K}
}
}
\references{
  High Performance Statistical Computing Website:
  \url{http://thirteen-01.stat.iastate.edu/snoweye/hpsc/}
}
\author{
  Wei-Chen Chen \email{wccsnow@gmail.com} and George Ostrouchov.
}
\seealso{
  \code{\link{em.step.worker}},
  \code{\link{aecm.step.worker}},
  \code{\link{apecm1.step.worker}}, \cr
  \code{\link{apecm2.step.worker}}.
}
\examples{
\dontrun{
# Examples can be found in the help pages of em.step.worker(),
# aecm.step.worker(), apecm1.step.worker(), and apecm2.step.worker().
}
}
\keyword{programming}

