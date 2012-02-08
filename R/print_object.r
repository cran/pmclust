### This files provides functions for output.

mb.print <- function(PARAM, CHECK){
  if(COMM.RANK == 0){
#    cat("=====\n")
    cat("\n")
    cat("Method: ", CHECK$method, "\n", sep = "")
    cat("Convergence: ", CHECK$convergence,
        "  iter: ", CHECK$iter,
        "  abs.err: ", CHECK$abs.err,
        "  rel.err: ", CHECK$rel.err, "\n", sep = "")
    cat("logL: ", PARAM$logL, "\n", sep = "")
    cat("K: ", PARAM$K, "\n", sep = "")
    if(CHECK$method %in% c("em", "aecm", "apecm1", "apecm2")){
      cat("\nETA:\n")
      print(PARAM$ETA)
    }
    cat("\nMU:\n")
    print(PARAM$MU)
    if(CHECK$method %in% c("em", "aecm", "apecm1", "apecm2")){
      cat("\nSIGMA:\n")
      print(matrix(do.call("c", PARAM$SIGMA), ncol = PARAM$K))
    }
#    cat("=====\n")
    cat("\n")
  }
  invisible(mpi.barrier())
} # End of mb.print().


### Print functions.
catmpi <- function(..., COMM.SHOW = 0){
  my.rank <- mpi.comm.rank()
  if(my.rank %in% COMM.SHOW){
    if(length(COMM.SHOW) > 1 || my.rank != 0){
      cat("== COMM.RANK: ", my.rank, "\n", sep = "")
    }
    cat(...)
  }
  invisible(mpi.barrier())
} # End of catmpi().

printmpi <- function(..., COMM.SHOW = 0){
  my.rank <- mpi.comm.rank()
  if(my.rank %in% COMM.SHOW){
    if(length(COMM.SHOW) > 1 || my.rank != 0){
      cat("== COMM.RANK: ", my.rank, "\n", sep = "")
    }
    print(...)
  }
  invisible(mpi.barrier())
} # End of printmpi().

