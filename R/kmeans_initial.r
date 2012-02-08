### This file gives a simple initialization.

initial.center.worker <- function(PARAM, MU = NULL){
  if(is.null(MU)){
    N.worker <- nrow(X.worker)
    N.allworkers <- mpi.allgather(as.integer(N.worker), type = 1,
                                  integer(COMM.SIZE))

    center.worker <- rep(0, PARAM$K)
    if(COMM.RANK == 0){
      center.worker <- sample(1:COMM.SIZE, PARAM$K, replace = TRUE,
                              prob = N.allworkers / PARAM$N) - 1
    }
    center.worker <- mpi.bcast(as.integer(center.worker), type = 1)

    tmp <- NULL
    n.center.worker <- sum(center.worker == COMM.RANK)
    if(n.center.worker > 0){
      id.center.worker <- sample(1:N.worker, n.center.worker)
      tmp <- matrix(X.worker[id.center.worker,], ncol = ncol(X.worker),
                    byrow = TRUE)
    }

    PARAM$MU <- unlist(mpi.allgather.Robj(obj = tmp))
    dim(PARAM$MU) <- c(PARAM$p, PARAM$K)
  } else{
    PARAM$MU <- MU
  }

  for(i.k in 1:PARAM$K){
    B <- W.plus.y(X.worker, -PARAM$MU[, i.k], nrow(X.worker), ncol(X.worker))
    Z.worker[, i.k] <<- -rowSums(B * B)
  }

  CLASS.worker <<- apply(Z.worker, 1, which.max)

  PARAM
} # End of initial.center.worker().

