#' Run ImageJ FFT filter macro from R
#' 
#' Passes an ImageJ FFT filter on files matching \code{BF.*.tif} in the target directory.
#' 
#' The modified images are saved to a new subdirectory with default name: "filtered/" (hardcoded in the macro).
#' 
#' Run \code{???imagej_macro_run} for more info.
#' 
#' @param pic.path Path to the directory containing the image files (passed as extra.args to imagej_macro_run).
#' @inheritParams imagej_macro_run
#' @inheritDotParams imagej_macro_run
#' @export
#' 
imagej_fft_filter <- function(
  pic.path,
  script.path = system.file("imagej_macros/FFT_filter_on_BFs_R.txt", 
                            package = "rcell2.cellid"),
  ...
  ){
  
  is.dir <- file.info(pic.path)$isdir
  
  if(is.dir) {
    pic.path <- paste0(normalizePath(pic.path), "/")
  } else {
    stop("The provided 'pic.path' is not a directory.")
  }
  
  imagej_macro_run(script.path = script.path, extra.args = pic.path, ...)
}

# #' Run ImageJ open on a path
# #' 
# #' @param pic.path Path to the image file.
# #' @inheritParams imagej_macro_run.headless
# #' @export
# #' 
# imagej.open <- function(
#   pic.path,
#   script.path = system.file("inst/imagej_macros/open.ijm", 
#                             package = "rcell2"),
#   wait = F){
#   
#   is.dir <- file.info(pic.path)$isdir
#   
#   if(!is.dir) {
#     pic.path <- paste0(normalizePath(pic.path))
#   } else {
#     stop("The provided 'pic.path' is a directory.")
#   }
#   
#   imagej_macro_run(script.path, pic.path)
# }


#' Run headless ImageJ Macro file
#' 
#' The macro program takes one argument: a path to a directory containing images matching the pattern "BF.*.out.tif".
#' 
#' Images are then proceesed with ImageJ's "Bandpass filter" and the output is saved to a subdirectory named "filtered".
#' 
#' @param imagej.path Path to the ImageJ binary (a path to "ImageJ-linux64" or equivalent).
#' @param script.path Path to the ImageJ macro. Defaults to built-in macro.
#' @param extra.args A string with extra arguments to the ImageJ command, pasted at the end.
#' @param headless Wether ImageJ should be run headlessly (no GUI).
#' @inheritParams base::system
#' 
imagej_macro_run <- function(
  script.path,
  imagej.path = "~/Software/FIJI/Fiji.app/ImageJ-linux64",
  wait = T, 
  headless = T,
  extra.args = ""){
  command <- paste(
    normalizePath(imagej.path),
    {if (headless) "--headless -macro" else "-macro"},
    normalizePath(script.path),
    extra.args
  )
  
  base::system(command, wait = wait)
}

