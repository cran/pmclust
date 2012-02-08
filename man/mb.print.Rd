\name{mb.print}
\alias{mb.print}
\title{Print results of model-based clustering}
\description{
  This function will print summarized messages for model-based clustering.
}
\usage{
  mb.print(PARAM, CHECK)
}
\arguments{
  \item{PARAM}{a set of convergent parameters to be printed.}
  \item{CHECK}{a set of checking parameters to be printed.}
}
\details{
  This function will provide a quick summary from the \code{PARAM} and
  \code{CHECK} typically the output of clusterings when algorithms stop.
  The \code{\link{COMM.RANK}} 0 will be used to take care the printing.
}
\value{
  Summarized messages will print/cat on screen by default.
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

