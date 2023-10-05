### Apply FFT filter on BF images

# Optionally, use ImageJ's bandpass FFT filter on the defocused brightfield images.
# 
# This can help Cell-ID find cells and reduce background artifacts.
# 
# List the BF images from the `arguments` function:

# Get renamed image paths
# data.dir <- file.path(data.dir, "renamed")
# 
# cellid.args <-
#   rcell2.cellid::arguments(path = data.dir,
#                            file.pattern = "^(BF|[TRY]FP)_Position(\\d+)_time(\\d+).tif$") |> 
#   filter(pos %in% 7:8)

#' Run Imagej's FFT filter on BF files
#' 
#' @param data.dir Data directory of the original images, appropriately named.
#' @param cellid.args Arguments dataframe as produced by \code{rcell2.cellid::arguments} at \code{data.dir}.
#' @param imagej.path Path to the ImageJ/FIJI binary.
#' @import foreach doParallel parallel
#' @export
run_fft_filter_on_bfs <- function(data.dir, 
                                  cellid.args, 
                                  imagej.path="~/Software/ImageJ/Fiji.app/ImageJ-linux64"){
  #### Split and symlink BFs
  # Split the BF files by position, and symlink them to different subdirectories:
  cellid.args.split <- split(cellid.args, ~pos)
  
  # Prefix or full path to where filtered images will be stored.
  fft.subdirs.prefix <- "bf_subset"
  
  fft.bfs.subdirs <- list()
  for(i in seq_along(cellid.args.split)) {
    # Save BF files
    bf.files <- cellid.args.split[[i]] %>%
      select(bf, path) %>% with(paste0(path,"/",bf)) %>% unique()
    
    # Make a directory to put them in
    fft.bfs.subdirectory <- paste0(fft.subdirs.prefix, "_", i)
    fft.bfs.dir <- file.path(data.dir, fft.bfs.subdirectory)
    dir.create(fft.bfs.dir)
    
    # Link them over
    file.symlink(from = normalizePath(bf.files),
                 to = file.path(dirname(bf.files[1]), fft.bfs.subdirectory, basename(bf.files)))
    
    fft.bfs.subdirs[i] <- fft.bfs.dir
  }
  
  # fft.bfs.subdirs
  
  #### Run filter
  
  # Run FFT filter in parallel, with one thread per subdirectory:
  
  n_cores <- parallel::detectCores() #- 1
  n_chunks <- length(fft.bfs.subdirs)
  
  # Make cluster
  cl <- parallel::makeCluster(
    spec = min(n_chunks,n_cores), 
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
  filtered.bfs.dir <- file.path(data.dir, "fft_images_dataset")
  dir.create(filtered.bfs.dir)
  
  # Symlink BF images to the final location.
  result <- file.symlink(from = normalizePath(filtered.bfs), 
                         to = file.path(filtered.bfs.dir, basename(filtered.bfs)))
  
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