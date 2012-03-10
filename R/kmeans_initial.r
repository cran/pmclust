### This file gives a simple initialization.

initial.center.spmd <- function(PARAM, MU = NULL){
  if(is.null(MU)){
    N.spmd <- nrow(X.spmd)
    N.allspmds <- mpi.allgather(as.integer(N.spmd), type = 1,
                                  integer(COMM.SIZE))

    center.spmd <- rep(0, PARAM$K)
    if(COMM.RANK == 0){
      center.spmd <- sample(1:COMM.SIZE, PARAM$K, replace = TRUE,
                              prob = N.allspmds / PARAM$N) - 1
    }
    center.spmd <- mpi.bcast(as.integer(center.spmd), type = 1)

    tmp <- NULL
    n.center.spmd <- sum(center.spmd == COMM.RANK)
    if(n.center.spmd > 0){
      id.center.spmd <- sample(1:N.spmd, n.center.spmd)
      tmp <- matrix(X.spmd[id.center.spmd,], ncol = ncol(X.spmd),
                    byrow = TRUE)
    }

    PARAM$MU <- unlist(mpi.allgather.Robj(obj = tmp))
    dim(PARAM$MU) <- c(PARAM$p, PARAM$K)
  } else{
    PARAM$MU <- MU
  }

  for(i.k in 1:PARAM$K){
    B <- W.plus.y(X.spmd, -PARAM$MU[, i.k], nrow(X.spmd), ncol(X.spmd))
    Z.spmd[, i.k] <<- -rowSums(B * B)
  }

  CLASS.spmd <<- apply(Z.spmd, 1, which.max)

  PARAM
} # End of initial.center.spmd().

