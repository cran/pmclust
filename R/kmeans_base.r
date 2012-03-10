### This file provides functions for kmeans.

kmeans.e.step.spmd <- function(PARAM){
  nrow <- nrow(X.spmd)
  ncol <- ncol(X.spmd)

  for(i.k in 1:PARAM$K){
    B <- W.plus.y(X.spmd, -PARAM$MU[, i.k], nrow, ncol)
    Z.spmd[, i.k] <<- sqrt(rowSums(B * B))
  }
} # End of kmeans.e.step.spmd().

kmeans.m.step.spmd <- function(PARAM){
  for(i.k in 1:PARAM$K){
    id <- CLASS.spmd == i.k
    tmp.n.id <- as.double(sum(id))
    tmp.n.id <- mpi.allreduce(tmp.n.id, type = 2, op = "sum")

    if(tmp.n.id > 0){
      tmp.sum <- colSums(matrix(X.spmd[id, ], ncol = PARAM$p))
    } else{
      tmp.sum <- rep(0.0, PARAM$p)
    }
    tmp.sum <- mpi.allreduce(tmp.sum, type = 2, op = "sum")

    PARAM$MU[, i.k] <- tmp.sum / tmp.n.id
  } 

  PARAM
} # End of kmeans.m.step.spmd().

kmeans.logL.step <- function(){
  tmp <- apply(Z.spmd, 1, which.min)
  tmp.diff <- sum(CLASS.spmd != tmp)

  CLASS.spmd <<- tmp
  mpi.allreduce(as.integer(tmp.diff), type = 1, op = "sum")
} # End of kmeans.logL.step().

check.kmeans.convergence <- function(PARAM.org, PARAM.new, i.iter){
    abs.err <- PARAM.new$logL
    rel.err <- abs.err / PARAM.new$N
    convergence <- 0

    if(i.iter > CONTROL$max.iter){
      convergence <- 2
    } else if(abs.err == 0 || rel.err < CONTROL$rel.err){
      convergence <- 1
    }

    if(CONTROL$debug > 1){
      catmpi("  check.em.convergence:",
             " abs: ", abs.err,
             ", rel: ", rel.err,
             ", conv: ", convergence, "\n", sep = "")
    }

    list(method = CHECK$method,
         iter = i.iter, abs.err = abs.err, rel.err = rel.err,
         convergence = convergence)
} # End of check.kmeans.convergence().

kmeans.step.spmd <- function(PARAM.org){
  CHECK <<- list(method = "kmeans", i.iter = 0, abs.err = Inf, rel.err = Inf,
                 convergence = 0)
  i.iter <- 1
  PARAM.org$logL <- PARAM.org$N

  ### For debugging.
  if((!is.null(CONTROL$save.log)) && CONTROL$save.log){
    if(! exists("SAVE.iter", envir = .GlobalEnv)){
      SAVE.param <<- NULL
      SAVE.iter <<- NULL
      CLASS.iter.org <<- unlist(apply(Z.spmd, 1, which.min))
    }
  }

  repeat{
    ### For debugging.
    if((!is.null(CONTROL$save.log)) && CONTROL$save.log){
      time.start <- proc.time()
    }

    PARAM.new <- kmeans.onestep.spmd(PARAM.org)

    CHECK <<- check.kmeans.convergence(PARAM.org, PARAM.new, i.iter)

    if(CHECK$convergence > 0){
      break
    }

    ### For debugging.
    if((!is.null(CONTROL$save.log)) && CONTROL$save.log){
      tmp.time <- proc.time() - time.start

      SAVE.param <<- c(SAVE.param, PARAM.new)
      CLASS.iter.new <- unlist(apply(Z.spmd, 1, which.min))
      tmp <- as.double(sum(CLASS.iter.new != CLASS.iter.org))
      tmp <- mpi.allreduce(tmp, type = 2, op = "sum")
      tmp.all <- c(tmp / PARAM$N, PARAM.new$logL,
                   PARAM.new$logL - PARAM.org$logL,
                   (PARAM.new$logL - PARAM.org$logL) / PARAM.org$logL)
      SAVE.iter <<- rbind(SAVE.iter, c(tmp, tmp.all, tmp.time))
      CLASS.iter.org <<- CLASS.iter.new
    }

    PARAM.org <- PARAM.new
    i.iter <- i.iter + 1
  }

  PARAM.new
} # End of kmeans.step.spmd().

kmeans.onestep.spmd <- function(PARAM){
#  if(COMM.RANK == 0){
#    Rprof(filename = "kmeans.Rprof", append = TRUE)
#  }

  PARAM <- kmeans.m.step.spmd(PARAM)
  kmeans.e.step.spmd(PARAM)

#  if(COMM.RANK == 0){
#    Rprof(NULL)
#  }

  PARAM$logL <- kmeans.logL.step()

  if(CONTROL$debug > 0){
    catmpi(">>kmeans.onestep: ", format(Sys.time(), "%H:%M:%S"),
           ", iter: ", CHECK$iter, ", logL: ",
                       sprintf("%-30d", PARAM$logL), "\n", sep = "")
    if(CONTROL$debug > 10){
      mb.print(PARAM, CHECK)
    }
  }

  PARAM
} # End of kmeans.onestep.spmd().


kmeans.update.class.spmd <- function(){
  CLASS.spmd <<- apply(Z.spmd, 1, which.min)
} # End of kmeans.update.class.spmd().

