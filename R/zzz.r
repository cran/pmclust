.First.lib <- function(lib, pkg)
{
  library.dynam("pmclust", pkg, lib)
} # End of .First.lib()

.Last.lib <- function(libpath)
{
  library.dynam.unload("pmclust", libpath)
} # End of .Last.lib()

