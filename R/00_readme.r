readme <- function(){
cat('
### Prespecified variables:
# (identical)
#   CONTROL: list[3]
#     - max.iter: integer[1], number of maximum iterations.
#     - abs.err: double[1], absolute tolerance.
#     - rel.err: double[1], relative tolerance.
#     - debug: integer[1], level of debug.
#     - RndEM.iter: integer[1], number of iterations of RndEM.
#     - exp.min: double[1], minimum exponent for base e.
#     - exp.max: double[1], maximum exponent for base e.
#   COMM.SIZE: integer[1], total workers.
#   p.times.logtwopi: p * log(2 * pi), for log likelihood.
# (different)
#   COMM.RANK: integer[1], rank in the communicator.
#   X.worker: double[N.worker, p], data.
#   ID.sample.worker: integer[], sampling id point to the large dataset.

### Global variables:
# (different)
#   Z.worker: double[N.worker, K], posterior probability.
#   Z.colSums: double[K], sum of posterior probability.
#   W.worker: double[N.worker, K], conditional log posterior probability.
#   W.worker.rowSums: double[N.worker], log density for each observations.
#   U.worker: double[N.worker, K], W.worker plus log eta.
#   CLASS.worker: double[N.worker], classification of observations.
#   CHECK: list[4], for output.
#     - method: string[1], "em", "aecm", "apecm1", "apecm2", or "kmeans".
#     - i.iter: integer[1], current iteration.
#     - abs.err: double[1], current absolute tolerance.
#     - rel.err: double[1], current relative tolerance.
#     - convergence: integer[1], status of convergence.
#                    * 0: keep running.
#                    * 1: converge successfully.
#                    * 2: run out of iterations.
#                    * 3: logL decreasing.

### Local variables: subjected to updated within EM iterations.
# (identical/bcast)
#   PARAM: list, contains all prameters changed over iterations.
#     - N: integer[1], total number of observations.
#     - p: integer[1], variable dimension. (p > 1)
#     - K: integer[1], number of clusters.
#     - ETA: double[K], mixing proportion.
#     - log.ETA: double[K], log of mixing proportion.
#     - MU: double[p, K], centers.
#     - SIGMA: list[K], dispersions.
#              - SIGMA[[k]]: double[p, p], for each component.
#                            * convert to double[p * (p + 1) / 2] for LAPACK.
#     - logL: double[1], current log likelihood.
#     - min.N.CLASS: integer[1], minimum size of cluster.

### Output: 
# (different)
#   CLASS.worker
# (rank = 0, only, identical)
#   PARAM, CHECK
')
} # End of readme().
