\name{pmclust-package}
\alias{pmclust-package}
\docType{package}
\title{Parallel Model-Based Clustering}
\description{
  The pmclust aims to utilize model-based clustering (unsupervised)
  for high dimensional and ultra large data, especially in a distributed
  manner. The package employs Rmpi to perform a parallel version of
  expectation and maximization (EM) algorithm for finite mixture
  Gaussian models. The unstructured dispersion matrices
  are assumed in the Gaussian models. The implementation is default in
  the single program multiple data (SPMD) programming model.
  The code can be executed through Rmpi and independent to most MPI
  applications. See the High Performance Statistical Computing (HPSC)
  website for more information, documents and examples.
}
\details{
\tabular{ll}{
  Package: \tab pmclust\cr
  Type: \tab Package\cr
%  Version: \tab 0.1-2\cr
%  Date: \tab 2010-02-18\cr
  License: \tab GPL\cr
  LazyLoad: \tab yes\cr
}

  The main function is \code{\link{em.step.spmd}} implementing the
  parallel EM algorithm for mixture multivariate Gaussian models
  with unstructured dispersions. This function groups a data matrix
  \code{\link{X.spmd}} into \code{K} clusters where \code{\link{X.spmd}} is
  potentially huge and taken from the global environment
  \code{\link[base]{.GlobalEnv}}.

  Other main functions \code{\link{aecm.step.spmd}},
  \code{\link{apecm1.step.spmd}}, and
  \code{\link{apecm2.step.spmd}} may provide better performance than the
  \code{\link{em.step.spmd}} in terms of computing time and convergent
  iterations.

  \code{\link{kmeans.step.spmd}} provides the fastest clustering among above
  algorithms, but it is restricted by Euclidean distance and spherical
  dispersions.
}
\references{
  High Performance Statistical Computing (HPSC) Website:
  \url{http://thirteen-01.stat.iastate.edu/snoweye/hpsc/}

  Chen, W.-C. and Maitra, R. (2011)
  \dQuote{Model-based clustering of regression time series data via
          APECM -- an AECM algorithm sung to an even faster beat},
  \emph{Statistical Analysis and Data Mining},
  \bold{4}, 567-578.

  Chen, W.-C. and Ostrouchov, G. (2012)
  \dQuote{Parallel Model-Based Clustering for Finite Mixture Gaussian Models},
  \emph{(in preparation)}.

  Dempster, A.P., Laird, N.M. and Rubin, D.B. (1977)
  \dQuote{Maximum Likelihood from Incomplete Data via the EM Algorithm},
  \emph{Journal of the Royal Statistical Society Series B},
  \bold{39}, 1-38.

  Lloyd., S. P. (1982)
  \dQuote{Least squares quantization in PCM},
  \emph{IEEE Transactions on Information Theory},
  \bold{28}, 129-137.

  Meng, X.-L. and Van Dyk, D. (1997)
  \dQuote{The EM Algorithm -- an Old Folk-song Sung to a Fast New Tune},
  \emph{Journal of the Royal Statistical Society Series B},
  \bold{59}, 511-567.
}
\author{
  Wei-Chen Chen \email{wccsnow@gmail.com} and George Ostrouchov
}
\seealso{
  \code{\link{em.step.spmd}},
  \code{\link{aecm.step.spmd}},
  \code{\link{apecm1.step.spmd}}, \cr
  \code{\link{apecm2.step.spmd}},
  \code{\link{kmeans.step.spmd}}.
}
\examples{
\dontrun{
# Under command mode, run the demo in 4 processors by
mpirun -np 4 Rscript -e \
  'demo(ex_em, package = "pmclust", ask = F, echo = F)'
mpirun -np 4 Rscript -e \
  'demo(ex_aecm, package = "pmclust", ask = F, echo = F)'
mpirun -np 4 Rscript -e \
  'demo(ex_apecm1, package = "pmclust", ask = F, echo = F)'
mpirun -np 4 Rscript -e \
  'demo(ex_apecm2, package = "pmclust", ask = F, echo = F)'
mpirun -np 4 Rscript -e \
  'demo(ex_kmeans, package = "pmclust", ask = F, echo = F)'
}
}
\keyword{package}

