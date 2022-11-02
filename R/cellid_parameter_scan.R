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
#' @param parameters.df Dataframe with one combination of parameters per row.
#' @param scan.arguments Output from \code{arguments}, filtered to your scanning needs.
#' @param test.dir Working directory for the parameter scan.
#' @param progress Print a progress bar. Requires the \code{doSNOW} package.
#' @inheritParams cell.load.alt
#' @inheritDotParams cell2
#'
#' @rawNamespace import(foreach, except = c("when", "accumulate"))
#' @return Data frame with the results. To use it, have a look at the rmarkdown template bundled in the package, or get it with \code{get_workflow_template_cellid()}.
#' @export
#'
#' @examples Have a look at the rmarkdown template bundled in the package, or get it with \code{get_workflow_template_cellid()}.
parameter_scan <- function(parameters.df, 
                           scan.arguments,
                           test.dir = normalizePath(paste0(tempdir(check = T), "/images_directory/test.dir"),
                                                    mustWork = F), 
                           fluorescence.pattern = "^([GCYRT]FP|[GCYRT]\\d+)_Position\\d+(?:_time\\d+)?.tif$",
                           progress = TRUE,
                           ...) {
  
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
                file.symlink(from = paste0(i["path"], "/",
                                           c(i["image"],i["bf"])),
                             to =   paste0(tmp.path, "/", 
                                           c(i["image"], i["bf"]))
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
              cellid.args.tmp$output <- paste(tmp.path, 
                                              basename(cellid.args.tmp$output),
                                              sep = "/")
              
              # Run Cell-ID
              cellid.out <- rcell2.cellid::cell2(arguments = cellid.args.tmp,
                                                 verbose = F, progress = F,
                                                 ...)
              
              # Load output
              cell.data <- rcell2.cellid::cell.load.alt(
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
        stack.name <- paste0(
          test.dir, "/", images$channel[1], "_stack", # "-pos_", images$pos[1], #"-time", images$t.frame[1],
          ".tif"
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
