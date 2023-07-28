.onLoad <- function(lib, pkg) {
  cat("Loading compiled code...\n")
  library.dynam("snglmm", pkg, lib)
}
