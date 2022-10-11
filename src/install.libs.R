execs <- c("cell")

if(.Platform$OS.type == "windows"){
  # https://stackoverflow.com/a/44273482/11524079
  execs <- paste0(execs, ".exe")
} 

if ( any(file.exists(execs)) ) {
  dest <- file.path(R_PACKAGE_DIR,  paste0('bin', R_ARCH))
  dir.create(dest, recursive = TRUE, showWarnings = FALSE)
  file.copy(execs, dest, overwrite = TRUE)
  unlink(execs)
}
