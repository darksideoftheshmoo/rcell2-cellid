#' Scan parameters for CellID
#' 
#' @details 
#' 
#' Requires the \code{doSNOW} package.
#' 
#' A few "summary" quantities are calculated for each run:
#' 
#' fft: mean fft.stat
#' 
#' el.p: mean el.p
#' 
#' ucids: Number of found cells
#' 
#' persistence: Average frames per cell
#' 
#' roughness: function copied from the seewave package: "Roughness or total curvature is the integrated squared second derivative" (see: https://cran.r-project.org/web/packages/seewave/index.html).
#' 
#' roughness2: This variation takes the sqrt of roughness and divides by total length. To correct for extra accumulated roughness in cells detected across more t.frames.
#'
#' Have a look at the examples in the rmarkdown template bundled in the package, or get it with \code{get_workflow_template_cellid()}.
#' 
#' @param parameters.df Dataframe with one combination of parameters per row.
#' @param scan.arguments Output from \code{arguments}, filtered to your scanning needs.
#' @param test.dir Working directory for the parameter scan. Creates a sub-directory of \code{tempdir()} if NULL (the default).
#' @param progress Print a progress bar if TRUE. Requires the \code{doSNOW} package.
#' @inheritParams get_cell_data
#' @inheritDotParams cell2
#'
#' @rawNamespace import(foreach, except = c("when", "accumulate"))
#' @return Data frame with the results. To use it, have a look at the rmarkdown template bundled in the package, or get it with \code{get_workflow_template_cellid()}.
#' @export
parameter_scan <- function(parameters.df, 
                           scan.arguments,
                           test.dir = NULL, 
                           fluorescence.pattern = "^([GCYRT]FP|[GCYRT]\\d+)_Position\\d+(?:_time\\d+)?.tif$",
                           progress = TRUE,
                           ...) {
  
  if(is.null(test.dir)){
    test.dir <- normalizePath(path = paste0(tempdir(check = T), "/images_directory/test.dir"),
                              mustWork = F)
  }
  
  # Record the amount of combinations
  test.params <- 1:nrow(parameters.df)
  test.pos <- unique(scan.arguments$pos)
  test.frames <- unique(scan.arguments$t.frame)
  
  # Clean test directory
  unlink(test.dir, recursive = T)
  dir.create(test.dir, recursive = T)
  
  ncores <- parallel::detectCores() - 1
  cl <- parallel::makeCluster(ncores, setup_strategy = "sequential")
  
  # Check if the progressbar can be shown
  if(!requireNamespace("doSNOW", quietly = T)){
    # If the doSNOW package is missing...
    if(progress) {
      # And if a progress-bar was requested...
      warning(paste0(
        "\nparameter_scan: a progress-bar was requested, but the doSNOW package is missing.",
        " Install the 'doSNOW' package to enable this feature.",
        " Now falling back to registerDoParallel, without progress-bar :(\n"
      ))
      # then fall back to no progress-bar.
      progress <- FALSE
    }
  }
  
  # Create a parallel back-end:
  if(progress){
    # Register doSNOW cluster
    doSNOW::registerDoSNOW(cl)
    # Setup a progressbar and run CellID:
    ntasks <- length(test.params)
    pb <- txtProgressBar(max = ntasks, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress=progress)
  } else {
    # Register a regular doParallel cluster
    doParallel::registerDoParallel(cl)
    opts <- list()
  }
  
  # Test
  test.param=test.params[1]
  
  # Run scan
  result <- 
    foreach(test.param=test.params, #.export = c(".parameters.df",)
            .options.snow=opts,
            .packages = c("base", "dplyr")) %dopar% {
              
              # Prepare a temp dir for each cellid run
              tmp.path <- tempfile(pattern = "dir", tmpdir = test.dir)
              dir.create(tmp.path)
              
              # Create shortcuts to the original images
              apply(scan.arguments, MARGIN = 1, function(i){
                # From the original path to the new temporary directory
                file.symlink(from = file.path(i["path"], c(i["image"], i["bf"])),
                             to =   file.path(tmp.path,  c(i["image"], i["bf"]))
                )
              })
              
              # Get one parameter set
              parameters.list.one <- parameters.df[test.param,]
              
              # Save parameters to tmpdir
              parameters.txt <- rcell2.cellid::parameters_write(parameters.list.one, 
                                                                param.dir = test.dir)
              
              # Copy original arguments
              cellid.args.tmp <- scan.arguments  # Assume it has been filtered already
              
              # Replace parameters with current ones
              cellid.args.tmp$parameters <- parameters.txt
              
              # Replace the original base paths
              cellid.args.tmp$path <- tmp.path
              cellid.args.tmp$output <- file.path(tmp.path, 
                                                  basename(cellid.args.tmp$output))
              
              # Run Cell-ID
              cellid.out <- rcell2.cellid::cell2(arguments = cellid.args.tmp,
                                                 verbose = F, progress = F,
                                                 ...)
              
              # Load output
              cell.data <- rcell2.cellid::get_cell_data(
                tmp.path,
                fluorescence.pattern = fluorescence.pattern
              )
              
              
              # Now compute some average "pseudo-metrics" about the segmentation.
              
              # "roughness" function copied from the seewave package
              # https://cran.r-project.org/web/packages/seewave/index.html
              roughness <- function(x, std=FALSE){
                if(std) x <- x/max(x)
                deriv2 <- diff(x, 1, 2)
                roughness <- sum(deriv2^2, na.rm=TRUE)
                return(roughness)
              } 
              # This variation takes the sqrt of roughness and divides by total length
              # To correct for extra accumulated roughness in cells detected
              # across more t.frames.
              roughness2 <- function(x, std=FALSE){
                if(std) x <- x/max(x)
                deriv2 <- diff(x, 1, 2)
                roughness <- sum(deriv2^2, na.rm=TRUE)
                roughness.norm <- sqrt(roughness)/length(x)
                # roughness <- sum(abs(deriv2), na.rm=TRUE)
                # roughness.norm <- (roughness)/length(x)
                return(roughness.norm)
              } 
              
              # Compute the summary numbers
              fft <- mean(cell.data$data$fft.stat)
              el.p <- mean(cell.data$data$el.p)
              ucids <- length(unique(cell.data$data$ucid)) # Number of found cells
              persistence <- nrow(cell.data$data) / ucids  # Average frames per cell
              roughness <- cell.data$data %>% group_by(ucid) %>% 
                summarise(roughness = roughness(a.tot)) %>% 
                with(roughness) %>% mean(na.rm = T)
              roughness2 <- cell.data$data %>% group_by(ucid) %>% 
                summarise(roughness2 = roughness2(a.tot)) %>% 
                with(roughness2) %>% mean(na.rm = T)
              
              
              # Keep only the image paths dataframe
              # And the "metrics"
              result <- 
                cell.data$images %>% filter(is.out) %>% 
                mutate(parameters = parameters.txt,
                       fft = fft,
                       el.p = el.p,
                       ucids = ucids,
                       persistence,
                       roughness,
                       roughness2
                )
              
              return(result)
            }
  
  # Close cluster
  close(pb)
  parallel::stopCluster(cl)
  
  
  # Bind the result:
  results.bound <- result %>% 
    bind_rows(.id = "id") %>%
    mutate(id = as.integer(id))
  
  # Prepare output  
  output_list <- list(
    results.bound=results.bound,
    test.params=test.params,
    test.pos=test.pos,
    test.frames=test.frames,
    test.dir=test.dir
  )
  
  # Chin-pum!
  return(output_list)
}


#' Make TIFF stacks for reviewing the result of parameter scans in ImageJ
#'
#' @param scan.results The full result from \code{parameter_scan}.
#' @param stack.channels Vector of channels which should be stacked.
#' @param annotation.font Font for the annotations. A mono-spaced font is recommended.
#'
#' @return A list with paths to the stacks, and 
#' @export
#' @importFrom stringr str_pad
#' @import dplyr
#'
make_scan_stacks <- function(scan.results, 
                             stack.channels = "BF.out", 
                             annotation.font = "Hack") {
  
  if(!requireNamespace("magick")){
    stop("make_scan_stacks: requires functions from the 'magick' package, which is not installed.")
  }
  
  # Load results
  test.dir <- scan.results$test.dir
  results.bound <- scan.results$results.bound
  test.params <- scan.results$test.params
  test.pos <- scan.results$test.pos
  test.frames <- scan.results$test.frames
  
  # Make stacks
  # images <- stack.paths[[1]]  # For testing
  stack.paths <- results.bound %>% 
    dplyr::filter(channel %in% stack.channels) %>% 
    dplyr::arrange(channel, id, t.frame, pos) %>% 
    # {split(., list(.$channel, .$pos))} %>%  # Fix for R's old split
    {split(., .$channel)} %>%  # Split only by channel
    lapply(function(images){
        
        # Prepare a file name for the stack
        stack.name <- file.path(test.dir, 
          paste0(
            images$channel[1], "_stack", # "-pos_", images$pos[1], #"-time", images$t.frame[1],
            ".tif"
          )
        )
        
        images.files <- images %>% 
          dplyr::arrange(channel, t.frame, pos, id) %>% 
          with(file)
        
        images.params <- images %>% 
          dplyr::arrange(channel, t.frame, pos, id) %>% 
          with(parameters) %>% 
          sapply(readr::read_file) %>% 
          unname()
        
        images.meta <- images %>% 
          dplyr::arrange(channel, t.frame, pos, id) %>% 
          with({
            paste0("channel: ", channel, "\n", 
                   "t.frame: ", t.frame, "\n", 
                   "position: ",     pos,     "\n")
          })
        
        images.info <- paste0(images.params, "\n", images.meta)
        
        images.info.padded <- sapply(images.info, function(images.info.pos){
          st <- strsplit(images.info.pos, "\\n")[[1]]
          mx <- max(nchar(st))
          st.padded <- stringr::str_pad(st, width = mx, side = "right")
          st.pasted <- paste(st.padded, collapse = "\n")
          return(st.pasted)
        })
        
        images.files %>% 
          magick::image_read() %>% 
          magick::image_annotate(images.info.padded, boxcolor = "white", color = "black", size = 10,
                                 font = annotation.font) %>% 
          magick::image_convert(colorspace = "Gray", matte = F) %>% 
          magick::image_write(path=stack.name, format = "tiff", depth = 8)
        
        return(stack.name)
      })
  
  ## Prepare ImageJ macros and print them to the console
  macros <- 
    lapply(stack.paths, function(stack.path){
      n_ch <- length(test.params)
      n_z <- length(test.pos)
      n_t <- length(test.frames)
      
      macro <- glue::glue(
        '// Macro',
        paste0('open("', stack.path, '");'),
        'run("Stack to Hyperstack...",',
        '"order=xyczt(default) channels=" + {n_ch} +',
        '" slices=" + {n_z} +',
        '" frames=" + {n_t} +',
        '" display=Grayscale");\n\n',
        .sep = "\n"
      )
    
    macro %>% cat()
    
    return(macro)
  })
  
  return(list(
    stack.paths=stack.paths,
    imagej.macros=macros
  ))
}



#' Prepare information text box to annotate an image with it's segmentation parameters and metadata
#' 
#' Internal use.
#' 
#' @param images The "results.bound" data frame produced by \code{parameter_scan}, filtered to contain only one channel.
#' @return A character vector for an information "text box", with each line padded to square length.
make_info_box <- function(images){
  
  if(length(unique(images$channel)) != 1) warning("make_info_box: multiple channel values were supplied.")
  
  # Generate contents of the text box with segmentation parameters.
  images.files <- images %>% 
    with(file)
  
  images.params <- images %>% 
    with(parameters) %>% 
    sapply(readr::read_file) %>% 
    unname()
  
  images.meta <- images %>% 
    with({
      paste0("channel: ", channel, "\n", 
             "t.frame: ", t.frame, "\n", 
             "position: ",     pos,     "\n")
    })
  
  images.info <- paste0(images.params, "\n", images.meta)
  
  images.info.padded <- sapply(images.info, function(images.info.pos){
    st <- strsplit(images.info.pos, "\\n")[[1]]
    mx <- max(nchar(st))
    st.padded <- stringr::str_pad(st, width = mx, side = "right")
    st.pasted <- paste(st.padded, collapse = "\n")
    return(st.pasted)
  })
  
  return(images.info.padded)
}


#' Annotate parameter scan output for reviewing the results in ImageJ
#' 
#' Warning: images in the scan result will be annotated and overwritten. Be sure to double check that your original "source" images are _not_ in the data frame (i.e. pass only output ".out.tif" images to this function).
#'
#' @param scan.results The full result from \code{parameter_scan}.
#' @param annotate.channels Vector of channels which should be annotated.
#' @param preserve_source_imgs If TRUE, an error will be raised when non-output (segmented) images are found in the input.
#' @param in.place If TRUE, the provided images will be annotated and overwritten. Else, the annotated images will be saved to a subdirectory of \code{test.dir} named by \code{annotated.imgs.dir}.
#' @param annotated.imgs.dir If \code{in.place} is FALSE, the annotated images will be saved to a sub-directory of \code{test.dir} named by \code{annotated.imgs.dir}, and suffixed by the channel names.
#' @param annotation.font Font for the annotations. A mono-spaced font is recommended.
#'
#' @return ImageJ macros to open the annotated images as a virtual stack.
#' @export
#' @importFrom stringr str_pad
#' @import dplyr
#'
annotate_scan_output <- function(scan.results, 
                                 annotate.channels = "BF.out", 
                                 preserve_source_imgs = TRUE,
                                 in.place = FALSE,
                                 annotated.imgs.dir = "annotated",
                                 annotation.font = "Hack") {
  
  if(!requireNamespace("magick")){
    stop("make_scan_stacks: requires functions from the 'magick' package, which is not installed.")
  }
  
  # Load results
  test.dir <- scan.results$test.dir
  results.bound <- scan.results$results.bound
  test.params <- scan.results$test.params
  test.pos <- scan.results$test.pos
  test.frames <- scan.results$test.frames
  
  # Check for source images in the input list.
  if(any(results.bound$is.out) & in.place & preserve_source_imgs) 
    stop("annotate_scan_output: non-output images detected for in place annotation. This would modify source data. To allow this, set 'preserve_non_out=FALSE'.")
  
  # Make stacks
  stack.paths.split <- results.bound %>% 
    dplyr::filter(channel %in% annotate.channels) %>% 
    dplyr::arrange(channel, id, t.frame, pos) %>% 
    # {split(., list(.$channel, .$pos))} %>%  # Fix for R's old split
    {split(., .$channel)}
  
  ## Annotate images, and prepare ImageJ macros and print them to the console
  macros <- stack.paths.split %>%  # Split only by channel
    lapply(function(images){
      # images <- stack.paths[[1]]  # For testing
      
      # Sort images.
      images.sorted <- images |> 
        dplyr::arrange(channel, t.frame, pos, id)
      # Get the file names and parameter set ID.
      images.files <- images.sorted %>% with(file)
      param.set.ids <- images.sorted %>% with(id)
      # Make info text boxes.
      images.info.padded <- make_info_box(images.sorted)
      
      # Make directory path for annotated images
      image.new.dir <- file.path(test.dir, paste0(annotated.imgs.dir, "-", images.sorted$channel[1])) |> 
        normalizePath(mustWork = F)
      
      # Read one image at a time, annotate it, and overwrite the original file.
      for (i in seq_along(images.files)) {
        image.file <- images.files[i]
        image.info.padded <- images.info.padded[i]
        param.set.id <- param.set.ids[i]
        
        image.annotated <- image.file %>% 
          magick::image_read() %>% 
          magick::image_annotate(text = image.info.padded, boxcolor = "white", color = "black", size = 10,
                                 font = annotation.font) %>% 
          magick::image_convert(colorspace = "Gray", matte = F)
        
        if(in.place){
          image.annotated %>% 
            magick::image_write(path=image.file, format = "tiff", depth = 8)
        } else {
          image.new.path <- file.path(image.new.dir, paste0(param.set.id, "-", basename(image.file)))
          dir.create(path = image.new.dir, showWarnings = FALSE)
          image.annotated %>% 
            magick::image_write(path=image.new.path, format = "tiff", depth = 8)
        }
      }
      
      # Prepare macro
      n_ch <- length(test.params)
      n_z <- length(test.pos)
      n_t <- length(test.frames)
      
      macro <- glue::glue(
        '// Macro',
        'run("Image Sequence...", "open={image.new.dir} file=(.*tif$) sort use");',
        'run("Stack to Hyperstack...", "order=xyczt(default) channels={n_t} slices={n_z} frames={n_ch} display=Grayscale"',
        ');\n\n',
        .sep = "\n"
      )
      
      macro %>% cat()
      
      return(macro)
    })
  
  return(macros)
}

