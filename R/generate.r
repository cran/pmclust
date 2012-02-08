### This file contains a simple data generation.

### N.K.worker: integer[K], number of elements of each cluster.
generate.basic.worker <- function(N.allworkers, N.worker, N.K.worker,
    N, p, K, seed){
  set.seed(seed)

  data.simu <- NULL
  data.class <- NULL
  data.n.class <- rep(0, K)
  for(i.k in 1:K){
    tmp.n.k <- N.K.worker[i.k]

    if(tmp.n.k > 0){
      tmp.data.simu <- NULL
      for(i.p in 1:p){
        mean <- i.k * 2 + i.p + 3
        sd <- sqrt(1 / mean)
        tmp.data.simu <- cbind(tmp.data.simu, rnorm(tmp.n.k, mean, sd))
      }
      data.simu <- rbind(data.simu, tmp.data.simu)
      data.class <- c(data.class, rep(i.k, tmp.n.k))
      data.n.class[i.k] <- tmp.n.k
    }
  }

  ret <- list(K = K, p = p, N = N, N.allworkers = N.allworkers,
              N.worker = N.worker, N.K.worker = N.K.worker, seed = seed,
              X.worker = data.simu, CLASS.worker = data.class,
              N.CLASS.worker = data.n.class)
  ret
} # End of generate.basic.worker().

