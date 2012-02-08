### This file gives initializations.

initial.em.worker <- function(PARAM, MU = NULL){
  if(is.null(MU)){
    ### Set semi-supervised information.
    N.worker <- nrow(X.worker)
    unlabeled.K <- PARAM$K
    unlabeled.N.worker <- N.worker
    # if(SS.clustering){
    #   unlabeled.K <- unlabeled.K - SS.K
    #   unlabeled.N.worker <- unlabeled.N.worker - length(SS.id.worker)
    # }

    ### Simple random sampling from data.
    N.allworkers <- mpi.allgather(as.integer(unlabeled.N.worker), type = 1,
                                  integer(COMM.SIZE))

    center.worker <- rep(0, unlabeled.K)
    if(COMM.RANK == 0){
      center.worker <- sample(1:COMM.SIZE, unlabeled.K, replace = TRUE,
                              prob = N.allworkers / sum(N.allworkers)) - 1
    }
    center.worker <- mpi.bcast(as.integer(center.worker), type = 1)

    tmp <- NULL
    n.center.worker <- sum(center.worker == COMM.RANK)
    if(n.center.worker > 0){
      N.pool <- 1:N.worker
      # if(SS.clustering && length(SS.id.worker) > 0){
      #   N.pool <- N.pool[-SS.id.worker]
      # }
      id.center.worker <- sample(N.pool, n.center.worker)
      tmp <- matrix(X.worker[id.center.worker,], ncol = ncol(X.worker),
                    byrow = TRUE)
    }
    MU <- matrix(unlist(mpi.allgather.Robj(obj = tmp)),
                 nrow = ncol(X.worker), ncol = unlabeled.K)

    ### Combind centers from semi-supervised information if any.
    # PARAM$MU <- cbind(SS.MU, MU)
    PARAM$MU <- MU
  } else{
    PARAM$MU <- MU
  }

  e.step.worker(PARAM)
  PARAM$logL <- logL.step()
  em.update.class.worker()

  PARAM
} # End of initial.em.worker().


initial.RndEM.worker <- function(PARAM){
  logL.save <- -Inf
  i.iter <- 1

  PARAM.org <- PARAM
  repeat{
    PARAM <- initial.em.worker(PARAM.org)
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
  PARAM <- initial.em.worker(PARAM.save, MU = PARAM.save$MU)
  PARAM
} # End of initial.RandRM.worker().

