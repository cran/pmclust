### This file gives initializations.

initial.em.spmd <- function(PARAM, MU = NULL){
  if(is.null(MU)){
    ### Set semi-supervised information.
    N.spmd <- nrow(X.spmd)
    unlabeled.K <- PARAM$K
    unlabeled.N.spmd <- N.spmd
    # if(SS.clustering){
    #   unlabeled.K <- unlabeled.K - SS.K
    #   unlabeled.N.spmd <- unlabeled.N.spmd - length(SS.id.spmd)
    # }

    ### Simple random sampling from data.
    N.allspmds <- mpi.allgather(as.integer(unlabeled.N.spmd), type = 1,
                                  integer(COMM.SIZE))

    center.spmd <- rep(0, unlabeled.K)
    if(COMM.RANK == 0){
      center.spmd <- sample(1:COMM.SIZE, unlabeled.K, replace = TRUE,
                              prob = N.allspmds / sum(N.allspmds)) - 1
    }
    center.spmd <- mpi.bcast(as.integer(center.spmd), type = 1)

    tmp <- NULL
    n.center.spmd <- sum(center.spmd == COMM.RANK)
    if(n.center.spmd > 0){
      N.pool <- 1:N.spmd
      # if(SS.clustering && length(SS.id.spmd) > 0){
      #   N.pool <- N.pool[-SS.id.spmd]
      # }
      id.center.spmd <- sample(N.pool, n.center.spmd)
      tmp <- matrix(X.spmd[id.center.spmd,], ncol = ncol(X.spmd),
                    byrow = TRUE)
    }
    MU <- matrix(unlist(mpi.allgather.Robj(obj = tmp)),
                 nrow = ncol(X.spmd), ncol = unlabeled.K)

    ### Combind centers from semi-supervised information if any.
    # PARAM$MU <- cbind(SS.MU, MU)
    PARAM$MU <- MU
  } else{
    PARAM$MU <- MU
  }

  e.step.spmd(PARAM)
  PARAM$logL <- logL.step()
  em.update.class.spmd()

  PARAM
} # End of initial.em.spmd().


initial.RndEM.spmd <- function(PARAM){
  logL.save <- -Inf
  i.iter <- 1

  PARAM.org <- PARAM
  repeat{
    PARAM <- initial.em.spmd(PARAM.org)
    if(class(PARAM) == "try-error"){
       next
    }

    N.CLASS <- get.N.CLASS(PARAM$K)
    if(any(N.CLASS < PARAM$min.N.CLASS)){
      next
    }

    if(CONTROL$debug > 0){
      catmpi("Initial: ", format(Sys.time(), "%H:%M:%S"),
             ", iter: ", i.iter, ", logL: ",
                         sprintf("%-20.10f", PARAM$logL), "\n", sep = "")
    }

    if(logL.save < PARAM$logL){
      logL.save <- PARAM$logL
      PARAM.save <- PARAM
      PARAM.save$initial.i.iter <- i.iter
    }

    i.iter <- i.iter + 1
    if(i.iter > CONTROL$RndEM.iter){
      break
    }
  }

  if(CONTROL$debug > 0){
    catmpi("Using initial iter: ", PARAM.save$initial.i.iter, "\n", sep = "")
  }
  PARAM <- initial.em.spmd(PARAM.save, MU = PARAM.save$MU)
  PARAM
} # End of initial.RandRM.spmd().

