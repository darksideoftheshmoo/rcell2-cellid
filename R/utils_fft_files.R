### Apply FFT filter on BF images

#' Run Imagej's FFT filter on BF files
#' 
#' Use ImageJ's bandpass FFT filter on the defocused brightfield images. This is optional, but can help Cell-ID find cells and reduce background artifacts.
#' 
#' Note that several intermediate files and directories are created during execution, which can remain there if the execution is interrupted for any reason.
#' 
#' To save disk space, fluorescence images are linked to the final directory as symbolic links instead of being duplicated by copying. Your milage may vary on a "non-unix OS" (i.e. Windows).
#' 
#' @param data.dir Data directory of the original images, appropriately named.
#' @param cellid.args Arguments dataframe as produced by \code{rcell2.cellid::arguments} at \code{data.dir}.
#' @param imagej.path Path to the ImageJ/FIJI binary.
#' @param fft.subdirs.prefix Name prefix for intermediate filtered BF files. 
#' @param filtered.bfs.dir.name Name for the final output directory,
#' @param n_cores Integer amout of processor cores to use.
#' @import foreach doParallel parallel
#' @export
run_fft_filter_on_bfs <- function(data.dir, 
                                  cellid.args, 
                                  fft.subdirs.prefix="bf_subset", 
                                  filtered.bfs.dir.name="fft_images_dataset",
                                  imagej.path="~/Software/ImageJ/Fiji.app/ImageJ-linux64",
                                  n_cores=NULL){
  
  # Default value for n_cores.
  if(is.null(n_cores)) n_cores <- max(1, parallel::detectCores() - 1)
  
  # Limit cores to the count of BFs if it is lower.
  n_cores <- length(unique(cellid.args$bf)) |> min(n_cores) |> max(1)
  
  #### Split and symlink BFs
  # Split the BF files by position, and symlink them to different subdirectories:
  # cellid.args.split <- split(cellid.args, ~pos)
  
  # Split the BF files by thread chunk, and symlink them to different subdirectories:
  cellid.args.split <- split(
    x = cellid.args, 
    f = rep(x = 1:n_cores, length.out = nrow(cellid.args))
  )
  
  if(!file.exists(imagej.path)) stop("Error! ImageJ executable not found at:" |> paste(imagej.path))
  
  # Prefix or full path to where filtered images will be stored.
  # fft.subdirs.prefix <- "bf_subset"
  
  fft.bfs.subdirs <- list()
  for(i in seq_along(cellid.args.split)) {
    # Save BF files
    bf.files <- cellid.args.split[[i]] %>%
      select(bf, path) %>% with(file.path(path, bf)) %>% unique()
    
    # Make a directory to put them in
    fft.bfs.subdirectory <- paste0(fft.subdirs.prefix, "_", i)
    fft.bfs.dir <- file.path(data.dir, fft.bfs.subdirectory)
    dir.create(fft.bfs.dir)
    
    # Link them over
    file.symlink(from = normalizePath(bf.files),
                 to = file.path(dirname(bf.files[1]), fft.bfs.subdirectory, basename(bf.files)))
    
    fft.bfs.subdirs[i] <- fft.bfs.dir
  }
  
  #### Run filter
  
  # Run FFT filter in parallel, with one thread per subdirectory:
  # n_chunks <- length(fft.bfs.subdirs)
  
  # Run FFT filter in parallel, with one subdirectory per thread:
  n_chunks <- n_cores
  
  # Make cluster
  cl <- parallel::makeCluster(
    spec = min(n_chunks, n_cores), 
    setup_strategy = "sequential"  # https://github.com/rstudio/rstudio/issues/6692
  )
  
  # This may take a while...
  # `%dopar%` <- foreach::`%dopar%`
  doParallel::registerDoParallel(cl)
  tryCatch(
    expr = {
      foreach(fft.bfs.subdir=fft.bfs.subdirs, .packages = "rcell2.cellid") %dopar% {
        imagej_fft_filter(pic.path = fft.bfs.subdir, 
                          imagej.path=imagej.path)
      }
    },
    error = function(e) {
      print("run_fft_filter_on_bfs error:")
      print(e)
      parallel::stopCluster(cl)
      stop("Couldn't run the ImageJ filter on one or more batches of Bf files.")
    })
  parallel::stopCluster(cl)
  
  #### Symlink filtered BFs
  
  # Symlink filtered BFs in each subdirectory to a single "filtered" images directory:
  
  filtered.bfs <- dir(file.path(fft.bfs.subdirs, "filtered"), full.names = T)
  
  # Set a name and make a path for the final location of the FFT-filtered images.
  filtered.bfs.dir <- file.path(data.dir, filtered.bfs.dir.name)
  dir.create(filtered.bfs.dir)
  
  # Move the BF images to the final location.
  result <- file.rename(from = normalizePath(filtered.bfs), 
                        # file.rename attempts to rename files (and from and to must be of the same length). 
                        # Where file permissions allow this will overwrite an existing element of to.
                        to = file.path(filtered.bfs.dir, basename(filtered.bfs)))
  
  # Remove intermediate directories.
  lapply(fft.bfs.subdirs, unlink, recursive=T) |> invisible()
  
  # all(result)  # Check
  
  # Symlink the rest of the images (fluorescence imgs) to the filtered directory:
  
  # Get all FL files (no BF files).
  fl.files <- cellid.args %>% 
    with(file.path(path, image)) %>% unique()
  
  # Symlink over to the directory of the FFT-filtered BFs
  result <- file.symlink(from = normalizePath(fl.files),
                         to = file.path(filtered.bfs.dir, basename(fl.files)))
  
  # all(result)  # Check
  
  # > Done!
  
  # You can run CellID using this directory as the `data.dir`:
  
  # "/mnt/arch_data/nicomic/Data/2023-08-18-Far1-Skot-240-to-3-nm/renamed/fft_images_dataset"
  new.data.dir <- normalizePath(filtered.bfs.dir)
  
  return(new.data.dir)
}