\name{catmpi}
\alias{catmpi}
\title{cat function}
\description{
  This function will cat messages from specified processors, by default
  messages will be shown on screen.
}
\usage{
  catmpi(..., COMM.SHOW = 0)
}
\arguments{
  \item{...}{a message to be printed.}
  \item{COMM.SHOW}{specified ranks to display messages.}
}
\details{
  \code{COMM.SHOW} can be a vector containing the ranks of processors which
  want to print messages.
}
\value{
  A \code{cat} will be called for the specified processors in \code{COMM.SHOW}
  and the messages will be shown on screen by default.
}
\references{
  High Performance Statistical Computing Website:
  \url{http://thirteen-01.stat.iastate.edu/snoweye/hpsc/}
}
\author{
  Wei-Chen Chen \email{wccsnow@gmail.com} and George Ostrouchov.
}
\seealso{
  \code{\link{printmpi}}.
}
\examples{
\dontrun{
### Save code in a file "demo.r" and run in 4 processors by
### > mpirun -np 4 Rscript demo.r

### Setup mpi environment.
library(Rmpi)
library(pmclust)
invisible(mpi.comm.dup(0, 1))

### Generate an example data.
COMM.RANK <- mpi.comm.rank()
catmpi("I am ", COMM.RANK, "\n", COMM.SHOW = c(1, 2))

### Quit Rmpi.
mpi.quit()
}
}
\keyword{programming}

