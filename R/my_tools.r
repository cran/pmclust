### This function initializes global variables.
set.global <- function(K = 2, PARAM = NULL,
    method = c("em", "aecm", "apecm1", "apecm2", "kmeans"),
    RndEM.iter = 10){
  ### Get data information.
  N.spmd <- nrow(X.spmd)
  N.allspmds <- mpi.allgather(as.integer(N.spmd), type = 1,
                                integer(mpi.comm.size()))
  N <- sum(N.allspmds)
  p <- ncol(X.spmd)

  ### Set parameters.
  if(is.null(PARAM)){
    PARAM <- list(N = N, p = p, K = K, ETA = rep(1 / K, K),
                  log.ETA = rep(-log(K), K), MU = NULL,
                  SIGMA = rep(list(diag(1.0, p)), K),
                  U = rep(list(), K),
                  U.check = rep(TRUE, K),
                  logL = NULL,
                  min.N.CLASS = (p + 1) * p * 0.5 + 1)
  } else{
    PARAM$N <- N
    K <- PARAM$K
  }

  ### Set global storages.
  CONTROL <<- list(max.iter = 1000, abs.err = 1e-4, rel.err = 1e-6,
                   debug = 1, RndEM.iter = RndEM.iter,
                   exp.min = log(.Machine$double.xmin),
                   exp.max = log(.Machine$double.xmax),
                   U.max = 1e+4,
                   U.min = 1e-6)

  COMM.SIZE <<- mpi.comm.size()
  COMM.RANK <<- mpi.comm.rank()

  p.times.logtwopi <<- p * log(2 * pi)

  Z.spmd <<- matrix(0.0, nrow = N.spmd, ncol = K)
  Z.colSums <<- rep(0.0, K)
  W.spmd <<- matrix(0.0, nrow = N.spmd, ncol = K)
  W.spmd.rowSums <<- rep(0.0, N.spmd)
  U.spmd <<- matrix(0.0, nrow = N.spmd, ncol = K)
  CLASS.spmd <<- rep(0, N.spmd)

  CHECK <<- list(method = method[1], i.iter = 0, abs.err = Inf, rel.err = Inf,
                 convergence = 0)

  ### For semi-supervised clustering.
#  assign.ss.spmd()

  for(i.k in 1:K){
    tmp.U <- decompsigma(PARAM$SIGMA[[i.k]])
    PARAM$U[[i.k]] <- tmp.U$value 
    PARAM$U.check[[i.k]] <- tmp.U$check 
  }

  PARAM
} # End of set.global().

