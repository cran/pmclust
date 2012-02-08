\name{assign.N.sample}
\alias{assign.N.sample}
\title{Obtain a set of random samples for X.worker}
\description{
  This utility function samples data randomly from \code{\link{X.worker}}
  to form a relatively small subset of original data. The EM algorithm on the
  smaller subset is topically performing fast and capturing rough structures of
  entire dataset.
}
\usage{
  assign.N.sample(total.sample = 5000, N.org.worker)
}
\arguments{
  \item{total.sample}{a total number of samples which will be selected from
                      the original data \code{\link{X.worker}}.}
  \item{N.org.worker}{the original data size,
                      i.e. \code{nrow(\link{X.worker})}.}
}
\details{
  This utility function performs simple random sampling without replacement
  for the original dataset \code{\link{X.worker}}.
  Different random seeds should
  be set before calling this function. A easy way is to call for example
  \code{\link[base]{set.seed}(123 + \link[Rmpi]{mpi.comm.rank}())}.
}
\value{
  A list variable will be returned and containing:
\tabular{ll}{\code{N}            \tab total sample size across all
                                      \eqn{S}{S} processors \cr
             \code{N.worker}     \tab sample size of given processor \cr
             \code{N.allworkers} \tab a collection of sample sizes for all
                                      \eqn{S}{S} processors \cr
             \code{ID.worker}    \tab index of selected samples ranged from 1
                                      to \code{N.org.worker}
}
  Note that \code{N} and \code{N.allworkers} are the same across all
  \eqn{S}{S} processors, but \code{N.worker} and \code{ID.worker} are most
  likely all distinct. The lengths of these elements are \eqn{1}{1} for
  \code{N} and \code{N.worker}, \eqn{S}{S} for \code{N.allworker}, and
  \code{N.worker} for \code{ID.worker}.
}
\references{
  High Performance Statistical Computing Website:
  \url{http://thirteen-01.stat.iastate.edu/snoweye/hpsc/}
}
\author{
  Wei-Chen Chen \email{wccsnow@gmail.com} and George Ostrouchov.
}
\seealso{
  \code{\link{set.global}}
}
\examples{
\dontrun{
# Save code in a file "demo.r" and run in 4 processors by
# > mpirun -np 4 Rscript demo.r

### Setup mpi environment.
library(Rmpi)
invisible(mpi.comm.dup(0, 1))

### Generate an example data.
set.seed(123 + mpi.comm.rank())
N.org.worker <- 5000 + sample(1:1000, 1)
ret.worker <- assign.N.sample(total.sample = 5000, N.org.worker)
cat("Rank:", mpi.comm.rank(), " Size:", ret.worker$N.worker,
    "\n", sep = "")

### Quit Rmpi.
mpi.quit()
}
}
\keyword{programming}

