### Controls

.PMC.CT <- list(
  algorithm = c("em", "aecm", "apecm", "apecma", "kmeans"),
  algorithm.gbd = c("em", "aecm", "apecm", "apecma", "kmeans"),
  method.own.X = c("gbdr", "spmdr", "common", "single"),
  CONTROL = list(
              max.iter = 1000,
              abs.err = 1e-4,
              rel.err = 1e-6,
              debug = 1,
              RndEM.iter = 10, 
              exp.min = log(.Machine$double.xmin),
              exp.max = log(.Machine$double.xmax),
              U.max = 1e+4,
              U.min = 1e-6,
              stop.at.fail = TRUE
            )
)

