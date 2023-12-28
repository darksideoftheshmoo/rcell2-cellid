#' Cell-ID output descriptions
#' 
#' @export
#' @param list.output Return the descriptions as a named list (TRUE), or as a data.frame (FALSE).
#' 
cellid_output_descriptions <- function(list.output=T){
  descs <- read.csv(sep = "\t", 
                    file = system.file("output_descriptions2.txt", package = "rcell2.cellid"))
  
  descs <- descs[!(is.na(descs$Variable.Name) | descs$Variable.Name == ""), -1:-2]
  
  if(list.output){
    descs.list <- as.list(descs$Description)
    names(descs.list) <- descs$Variable.Name
  return(descs.list)
  } else {
    return(descs)
  }
}

#' Cell-ID output description
#' @export
output_help <- cellid_output_descriptions(list.output = T) |> suppressWarnings()


#' Cell-ID output description
#' @export
output_help_df <- cellid_output_descriptions(list.output = F) |> suppressWarnings()

#' Cell-ID parameter descriptions
#' 
#' @param list_format If TRUE then format the dataframe into a named list
#' @export
#' 
cellid_parameter_descriptions <- function(list_format=T){
  
  warning("Warning: not all Cell-ID input parameters are fully documented. Refer to the original publication for further detail.")
  
  descs <- read.csv(sep = "\t", 
                    file = system.file("parameters_description.tsv", package = "rcell2.cellid"))
  
  descs <- descs[!(is.na(descs[[1]]) | descs[[1]] == ""),]
  
  if(!isTRUE(list_format)) return(descs)
  
  descs.split <- split(descs, ~parameter)
  descs.list <- lapply(descs.split, function(d) setNames(c(d), names(d)))
  
  return(descs.list)
}

#' Cell-ID parameter descriptions
#' @export
parameters_help <- cellid_parameter_descriptions(list_format = T) |> suppressWarnings()


#' Cell-ID parameter descriptions
#' @export
parameters_help_df <- cellid_parameter_descriptions(list_format = F) |> suppressWarnings()


# #' Correr Cell-ID desde R usando .C()
# #'
# #' @param args el comando de cellid entero, tal como se ejecutaria en bash "cell -p ..."
# #' @param debug_flag Set to 0 to disable CellID printf messages.
# # @useDynLib rcell2 CellID
# #' @export
# #' @return Exit code from CellID
# cellid <- function(args, debug_flag=0){
# 
#   # args <- "~/Software/cellID-linux/cell -p ~/Projects/Colman/HD/scripts/cellMagick/data/images/parameters.txt -b /tmp/Rtmp7fjlFo/file2b401093d715 -f /tmp/Rtmp7fjlFo/file2b402742f6ef -o ~/Projects/Colman/HD/uscope/20200130_Nico_screen_act1_yfp/1/Position001/out"
#   argv <- strsplit(args, " ")[[1]]
#   argc <- length(argv)
#   
#   if(debug_flag != 0) print("Printing argv and argc before .C() call to CellID.")
#   if(debug_flag != 0) print(argv)
#   if(debug_flag != 0) print(argc)
# 
#   exit_code <- 0
#   
#   if(exit_code != 1) stop(paste("CellID is not bundled in this branch, see master_cellid", exit_code))
#   
#   return(exit_code)
# }

#' Path to the installed cell binary
#' 
cell2_command <- function(){
  system.file(c("bin/cell", "bin/x64/cell.exe"), package = "rcell2.cellid", mustWork = T)[1]
}

#' Test the installed cell binary
#' @inheritDotParams system
cell2_test <- function(...){
  system(command = paste(cell2_command(), "-h"), ...)
}

#' Function to run Cell-ID
#' 
#' Handles Cell-ID executions: creates necessary directories, sets-up parallelization, wraps the Cell-ID binary, and checks exit status.
#'
#' @param arguments An argument data.frame, as built by \code{rcell2.cellid::arguments}.
#' @param cell.command By default \code{NULL}, to use the built-in binary. Otherwise a path to a Cell-ID binary executable (get if from https://github.com/darksideoftheshmoo/cellID-linux).
#' @param n_cores Number of cores to use for position-wise parallelization,internally capped to number of positions in \code{arguments}. Set to 1 to disable parallelization. If NULL, defaults to available cores - 1.
#' @param dry Do everything without actually running Cell-ID, print the commands that would have been issued.
#' @param debug_flag Set to 0 to disable Cell-ID printf messages (built-in Cell-ID only).
#'
#' @param output_coords_to_tsv Set to TRUE to write cell interior and boundary pixels coordinates of each cell to a compressed \code{.tsv} file, located in the main output directory (Cell-ID option '-t'). This data can be loaded with the \code{cell.load.boundaries} function.
#'
#' @param encode_cellID_in_pixels Set to TRUE to encode cellIDs in the intensity values of the boundary pixels, and blank the rest of the output image (Cell-ID option '-m'). Pixel intensities are proportional to each cellID, following the relationship \code{cellID = 65535 - boundary_intensity - 1}. Only boundary pixels are used by default; this behavior can be modified by enabling \code{label_cells_in_bf}, \code{fill_interior_pixels}, or \code{interior_offset}.
#'
#' @param label_cells_in_bf Set mask boundary pixel intensities proportional to each cellID, and add cellID numbers to the cells with maximum pixel intensity \code{65535} (Cell-ID option '-l', default FALSE).
#' 
#' @param fill_interior_pixels Fill each cell interior area in the output BF image file with intensity-labeled pixels (Cell-ID option '-i'). This overrides cell labeling.
#'
#' @param interior_offset Offset boundary and interior pixel intensities by a calculated 'interior_offset' threshold. cellID will relate to interior pixel intensities with the relationship \code{cellID = 65535 - boundary_intensity - interior_offset - 1}. The offset defaults to 5000, but may have a larger value for images or time series with more than 2500 cells (Cell-ID option '-w').
#'
#'
#' @param write_initial_time Write the absolute time of the first image to a text file (Cell-ID option '-z').
#' @param save.logs Set to TRUE to save Cell-ID logs to text files, into the output directory of their corresponding position.
#' @param verbose Print start-up messages.
#' @param progress Print a progress bar. Requires the \code{doSNOW} package.
#' @inheritParams arguments
#' 
#' @return A data.frame with one column indicating the issued commands and exit codes (in the command.output column). If the execution was successful, now you may run \code{rcell2::load_cell_data} or \code{rcell2.cellid::get_cell_data} to get the results from the Cell-ID output, typically located at the images path.
#' 
#' @import purrr dplyr stringr tidyr doParallel readr parallel
#' @importFrom foreach foreach %dopar% %do%
#' @rawNamespace import(foreach, except = c("when", "accumulate"))
#' @importFrom purrr map
#' @export
cell2 <- function(arguments,
                  cell.command = NULL,
                  n_cores = NULL, 
                  debug_flag=0,
                  dry = F,
                  output_coords_to_tsv = F,    # -t flag
                  encode_cellID_in_pixels = F, # -m flag
                  label_cells_in_bf = F,       # -l flag
                  fill_interior_pixels = F,    # -i flag
                  interior_offset = F,         # -w flag
                  write_initial_time = F,      # -z flag
                  save.logs = T, verbose=F,
                  progress=T,
                  check_fail=F){
  
  if(F){
    # For testing (NOT RUN)
    n_cores = 2
    debug_flag=0
    dry = F
    label_cells_in_bf = F
    fill_interior_pixels = F
    output_coords_to_tsv = F
    encode_cellID_in_pixels = F
    save.logs = T
  }
  
  # Check for existing output files
  arguments_check(arguments, check_fail)
  
  # Cell-ID path setup ####
  if(is.null(cell.command)){
    tryCatch(
      expr = {
        cat("\nUsing built-in Cell-ID binary.")
        if(verbose) cat(" Printing Cell-ID's help message:\n") else cat("\n")
        cell.command <- cell2_command()
        system(paste(cell.command, "-h"), ignore.stdout = verbose)
        if(verbose) cat("\n\n")
      },
      error = function(e) {
        print("cell2 error:")
        print(e)
        stop("Couldn't use the builtin Cell-ID, please specify a path to an external executable.")
      })
  } else {
    cat("\nUsing custom Cell-ID binary.")
    if(verbose) cat(" Printing Cell-ID's help message:\n") else cat("\n")
    warning("Custom Cell-ID binary selected. Keep in mind that rcell2.cellid::cell2 is meant to wrap the updated version. Double check the outputs.")
    system(paste(cell.command, "-h"), ignore.stdout = verbose)
    if(verbose) cat("\n\n")
  }
  # Check the path
  stopifnot(file.exists(cell.command))
  
  # Get some position and t.frame info
  positions <- arguments$pos %>% unique()
  n_positions <- positions %>% length()
  
  # Create output directories ####
  for(d in unique(arguments$output)) dir.create(d, showWarnings = F)
  
  # Prepare parallel backend if requested ####
  if(is.null(n_cores)) n_cores <- parallel::detectCores() - 1
  # Use parallel backend only if there are two or more cores
  if(n_cores >= 2){
    # Choose parallel foreach operator
    `%do_op%` <- foreach::`%dopar%`
    
    # Prepare cluster logfile
    cl.outfile <- tempfile(pattern = "dopar", 
                           # tmpdir = "/tmp", 
                           fileext = ".log")
    cat(paste("---- cell2: writing logs to file:", cl.outfile))
    
    # Make cluster
    cl <- parallel::makeCluster(
      spec = min(n_positions,n_cores), 
      outfile = cl.outfile, # outfile = NULL
      setup_strategy = "sequential"  #https://github.com/rstudio/rstudio/issues/6692
    )
    
    # Export arguments dataframe to cluster
    parallel::clusterExport(cl, "arguments", envir = environment())
    
    # Check progressbar config.
    if(progress & !requireNamespace("doSNOW", quietly = T)){
      warning("cell2: a progressbar was requrested but the 'doSNOW' is not installed. Disabling progress bar.")
      progress <- FALSE
    }
    
    # Register cluster
    if(!progress){
      doParallel::registerDoParallel(cl)
      opts <- list()
    } else {
      # Use registerDoSNOW for a progress bar
      doSNOW::registerDoSNOW(cl)
      # Setup the progressbar
      ntasks <- length(positions)
      pb <- txtProgressBar(max = ntasks, style = 3)
      progress <- function(n) setTxtProgressBar(pb, n)
      opts <- list(progress=progress)
    }
  } else {
    # Choose sequential foreach operator
    `%do_op%` <- foreach::`%do%`
  }
  
  # Run CellID ####
  # Iterate over positions
  sent_commands <- foreach::foreach(pos=positions, .options.snow=opts) %do_op% {
    # Get arguments for current position
    arguments_pos <- arguments[arguments$pos == pos,]
    
    # Get the path to the parameters
    parameters <- arguments_pos$parameters[1]
    
    # Make output path with file prefix
    output_prefix <- paste0(normalizePath(arguments_pos$output[1]), "/out")
    
    # Prepare text files paths
    bf_rcell2 <- tempfile(tmpdir = arguments_pos$output[1],
                          fileext = ".txt",
                          pattern = "bf_rcell2-")
    fl_rcell2 <- tempfile(tmpdir = arguments_pos$output[1],
                          fileext = ".txt",
                          pattern = "fl_rcell2-")
    dark <- tempfile(tmpdir = arguments_pos$output[1],
                     fileext = ".txt",
                     pattern = "dark-")
    flat <- tempfile(tmpdir = arguments_pos$output[1],
                     fileext = ".txt",
                     pattern = "flat-")
    
    # Image path vectors
    bf.imgs <- file.path(arguments_pos$path, arguments_pos$bf)
    fl.imgs <- file.path(arguments_pos$path, arguments_pos$image)
    
    # Check file existence
    stopifnot(all(
      file.exists(bf.imgs, fl.imgs)
    ))
    
    # Write image lists to text files
    base::write(x = bf.imgs, file = bf_rcell2)
    base::write(x = fl.imgs, file = fl_rcell2)
    
    # Write dark correction image lists to a text file
    use_dark <- "dark" %in% names(arguments_pos)
    if(use_dark) {
      if(verbose) cat(paste0("\nUsing dark image corrections from file: ", dark, "\n"))
      dark.imgs <- file.path(arguments_pos$path, arguments_pos$dark)
      base::write(x = dark.imgs, file = dark)
    }
    
    # Write flat correction image lists to a text file
    use_flat <- "flat" %in% names(arguments_pos)
    if(use_flat) {
      if(verbose) cat(paste0("\nUsing flat image corrections from file: ", flat, "\n"))
      flat.imgs <- file.path(arguments_pos$path, arguments_pos$flat)
      base::write(x = flat.imgs, file = flat)
    }
    
    # Third image support (nuclear/vacuole tagging)
    use_third <- "third" %in% names(arguments_pos)
    if(use_third){
      # TODO: force error if parameter not in config.
      third_rcell2 <- tempfile(tmpdir = arguments_pos$output[1],
                               fileext = ".txt",
                               pattern = "3rd_rcell2-")
      third.imgs <- file.path(arguments_pos$path, arguments_pos$third)
      base::write(x = third.imgs, file = third_rcell2)
    }
    
    # Check cell path again (not really needed i guess)
    if(is.null(cell.command)) stop("\nError: cell.command must point to an existing CellID binary on your system.")
    
    # Build the arguments for the Cell-ID command
    command.args <- paste0(
      # Mandatory arguments
      " -b ", bf_rcell2,
      " -f ", fl_rcell2,
      " -o ", output_prefix,
      {if(use_third) paste0(" -g ", third_rcell2) else ""},
      " -p ", parameters,
      # Correction files
      {if(use_dark) paste0(" -D ", dark) else ""},
      {if(use_flat) paste0(" -F ", flat) else ""},
      # Flags
      {if(label_cells_in_bf)       " -l" else ""},
      {if(output_coords_to_tsv)    " -t" else ""},
      {if(encode_cellID_in_pixels) " -m" else ""},
      {if(fill_interior_pixels)    " -i" else ""},
      {if(interior_offset)         " -w" else ""},
      {if(write_initial_time)      " -z" else ""}
    )
    
    # Paste to the binary's path
    command <- paste0(normalizePath(cell.command),
                      command.args
    )
    
    # Default is to print messages to console
    cellid.log <- ""
    cellid.err <- ""
    
    # Otherwise, write command and outputs to log files
    if(save.logs){
      # Save parameters
      cellid.pars <- arguments_pos$parameters[1]
      file.copy(from = cellid.pars, to = paste0(arguments_pos$output[1], .Platform$file.sep),
                overwrite = T)
      # Save cellid standard output
      cellid.log <- tempfile(tmpdir = arguments_pos$output[1],
                             fileext = ".txt",
                             pattern = "cellid_log-")
      # Save cellid standard error
      cellid.err <- tempfile(tmpdir = arguments_pos$output[1],
                             fileext = ".txt",
                             pattern = "cellid_error-")
      # Save cellid system command
      cellid.cmd <- tempfile(tmpdir = arguments_pos$output[1],
                             fileext = ".txt",
                             pattern = "cellid_cmd-")
      write(x = c("# Cell-ID command:\n\n", command, "\n\n"),
            file = cellid.cmd)
    } 
    
    # Skip execution if "dry" mode was selected
    if(!dry) {
      # Otherwise, execute CellID
      command.output <- system2(command = normalizePath(cell.command),
                                args = command.args,
                                stdout = cellid.log,
                                stderr = cellid.err,
                                wait = T)
      
    } else {
      command.output = NULL
    }
    
    cat(paste0("---- cell2: Done with position ", pos, "\n "))
    
    return(
      list(
        command.output = command.output,
        command = command,
        cellid.log = cellid.log,
        cellid.err = cellid.err
      )
    )
  }
  
  # Close cluster ####
  if(n_cores > 1){
    parallel::stopCluster(cl)
  }
  
  # Prepare output ####
  output <- dplyr::bind_rows(sent_commands)
  
  # Check if CellID's exit codes are weird
  if(!dry){
    if(all(output$command.output == 0)){
      cat("\nDone, please examine logs if anything seems strange :)")
    } else {
      cat("\nDone, please examine logs immediately! some exit codes signaled errors :(")
    }
  }
  
  return(output)
}

#' foreach and parLapply cluster test
#' 
#' @keywords internal
#' @import parallel doParallel
#' @importFrom foreach %dopar%
cluster_test <- function(){
  cl <- parallel::makeCluster(2)
  
  doParallel::registerDoParallel(cl)
  
  # Prueba con base
  parallel::parLapply(cl, list(1,2), function(x) print(x))
  
  # Prueba con foreach
  foreach::foreach(x=list(1,2)) %dopar% print(x)
  
  parallel::stopCluster(cl)
}

#' Obtener argumentos para CellID
#' 
#' @details 
#' 
#' All 4 regex groups are mandatory, but they may be left empty to exclude a field.
#' For instance, \code{t.frame} may be left as empty parenthesis \code{()},
#' but always preserving group order defined by \code{file.pattern.groups.order}.
#' 
#' If you only have BF images, consider setting \code{bf_as_fl = T}.
#' 
#' The "pos" and "channel" regex groups must always match position and channel identifiers in the file name, respectively.
#' 
#' Here are some example \code{file.pattern} regular expressions, when \code{file.pattern.groups.order = c("ch", "pos", "t.frame")}:
#' 
#' - No Z planes, no time (note the empty parentheses): \code{file.pattern = "^(BF|[TYR]FP)_Position(\\d+)().tif$"}
#' 
#' - No Z planes, with time: \code{file.pattern = "^(BF|[A-Z]FP)_Position(\\d+)_time(\\d+).tif$"}
#' 
#' - With Z planes and time: \code{file.pattern = "^(BF|[TYR]FP|[TYR]\\d{2})_Position(\\d+)_time(\\d+).tif$"}
#' 
#' Note that "Z planes" will be processed by CellID the same way as other channels.
#' More importantly, CellID distinguishes and groups image channels using the first 3 letters in the file name,
#' which limits planes to 100 per channel letter (R00 to R99).
#'
#' @param path Directory where images are stored, full path.
#' @param parameters Path to the parameters file or a data.frame with "pos" (position number) and "parameter" (path) columns. Defaults to \code{parameters_write()}.
#' @param BF.pattern Regex pattern to detect BF images only. Defaults to: \code{"^BF"}
#' @param file.pattern Regex pattern for all tif files, with one group for each of \code{c("ch", "pos", "t.frame")} in \code{file.pattern.groups.order}. Uses \code{"^(BF|[A-Z]FP)_Position(\\d+)_time(\\d+).tif$"} by default. To omit time, use an empty group for the t.frame in the regex, for example: \code{"^(BF|[A-Z]FP)_Position(\\d+)().tif$"}.
#' To consider Z-stacks, use something like "^(BF|[A-Z]\\d+)_Position(\\d+)_time(\\d+).tif$"
#' @param file.pattern.groups.order A character vector of components \code{c("ch", "z", "pos", "t.frame")} with order corresponding to the order of groups in \code{file.pattern}.
#' @param output.dir.basename Basename for the CellID output directories for each position.
#' @param tiff.ext Regex pattern for the tif file extension.
#' @param bf_as_fl If TRUE, BF paths will be used as FL paths. This allows for BF-only experiments.
#' @param dark_pattern A regular expression which matches only "dark correction" tiff files in the path; and has a single capturing group for the channel ID. If the exposure times match, Cell-ID will subtract this background image to fluorescent images (disregarding channel). Note: if the exposure of an FL image does not match any dark images, the correction is skipped without a warning.
#' @param flat_pattern A regular expression which matches only "flat correction" tiff files in the path; and has a single capturing group for the channel ID. With thses images, Cell-ID will attempt to correct for uneven illumination, on the corresponding fluorescent channels. Note: Cell-ID matches images to corrections using the three-letter prefix, and uses only the correction closest in time to the current FL image.
#' @inheritParams arguments_check
#' @return A data.frame with all the information needed to run CellID
#' @import dplyr tidyr
# @examples
# cell.args <- cellArgs(path = path)
#' @export
arguments <- function(path,
                      parameters=rcell2.cellid::parameters_write(),
                      BF.pattern = "^BF",
                      # file.pattern = "^(BF|[A-Z]FP)()_Position(\\d+)_time(\\d+).tif$",
                      # file.pattern.groups.order = c("ch", "z", "pos", "t.frame"),
                      file.pattern = "^(BF|[A-Z]FP)_Position(\\d+)_time(\\d+).tif$",
                      file.pattern.groups.order = c("ch", "pos", "t.frame"),
                      output.dir.basename = "Position",
                      # out.dir = "out",
                      tiff.ext = "tif$",
                      bf_as_fl=F,
                      dark_pattern="^(BF|[A-Z]FP)_dark.tif$",
                      flat_pattern="^(BF|[A-Z]FP)_flat.tif$",
                      check_fail = F
                      ){
  
  # Check capturing groups order argument
  if(!identical(sort(file.pattern.groups.order),
                sort(c("ch", "pos", "t.frame")))) 
    stop('arguments error: file.pattern.groups.order must contain "ch", "pos", and "t.frame" (in an appropriate order).')
  if (file.pattern.groups.order[1] != "ch") {
    stop('arguments error: the first item in file.pattern.groups.order must be "ch". Cell-ID uses the first three letters to group imaging channels.')
  }
  
  # Normalize the images path
  path <- normalizePath(path)
  
  # Make images dataframe ####
  # List images at the path
  pic_files <- dir(path, pattern = file.pattern)
  # Make a simple first check
  if(length(pic_files) == 0) stop(paste("arguments error: no image files retrieved using file pattern:", file.pattern))
  
  # Extract channle, position, and frame information.
  pics_df <- data.frame(image = pic_files,
                        path = path) %>% 
    tidyr::extract(col = image, 
                   into = file.pattern.groups.order,
                   regex = file.pattern, 
                   remove = F)
  
  # Get BFs: All BFs are BFs
  brihtfield_pics <- pics_df %>% 
    dplyr::filter(str_detect(string = image,
                             pattern = BF.pattern)) %>% 
    dplyr::rename(bf = image) %>% 
    dplyr::select(pos, t.frame, bf)
  
  # Get FLs: Use BFs as brightfields?
  if(isTRUE(bf_as_fl)){
    # Simple:
    fluor_pics <- pics_df
    
    # Warn user of potential unintended usage of option
    found_fl <- any(with(pics_df, str_detect(string = image,
                                             pattern = BF.pattern, 
                                             negate = T)))
    if(found_fl) warning("arguments warning: BF as FL option enabled, but some FL files were found. Please check if this is what you need.")
  } else {
    # Else, all non-BFs are fluor images:
    fluor_pics <- pics_df %>% 
      filter(str_detect(string = image,
                        pattern = BF.pattern, 
                        negate = T))
  }
  
  # Checks
  if(nrow(brihtfield_pics) == 0) stop("arguments error: Brightfield images missing, Check your directories and file.pattern.")
  if(nrow(fluor_pics) == 0) stop("arguments error: Fluorescence images missing, but BFs were found. Check your directories, file.pattern, and consider setting bf_as_fl.")
  
  # Bind df's
  arguments.df <- dplyr::left_join(
    fluor_pics,
    brihtfield_pics,
    by = c("pos", "t.frame")
  )
  
  # Check for missing BFs
  if(any(is.na(arguments.df$bf))){
    filter(arguments.df, is.na(bf)) %>% print()
    stop("arguments error: there are missing brightfield images")
  }
  
  # Add output column ####
  # and arrange by position and t.frame
  arguments.df.out <- arguments.df %>% 
    mutate(output = paste0(path, .Platform$file.sep, output.dir.basename, pos)) %>% 
    mutate(pos = as.integer(pos),
           t.frame = as.integer(t.frame)) %>% 
    arrange(pos, t.frame)
  
  # Add parameters column ####
  # Recycle parameters if lenght is 1
  if(length(parameters) == 1 & is.atomic(parameters)){
    arguments.df.out$parameters <- parameters
  } else {
  # Else bind to the passed parameters data.frame
    arguments.df.out <- left_join(arguments.df.out,
                                  dplyr::select(parameters, pos, parameters),
                                  by = "pos")
  }
  # Normalize parameters' paths
  arguments.df.out <- arguments.df.out %>% 
    mutate(parameters = normalizePath(parameters))
  
  # Get correction files ####
  
  # List dark files
  dark_files <- dir(path, pattern = dark_pattern)
  if(length(dark_files) > 0){
    # Prepare DF
    dark_channels <- data.frame(dark = dark_files,
                                ch = sub(dark_pattern, "\\1", dark_files))
    # Compare correction channels with image channels 
    if(!setequal(dark_channels$ch, arguments.df.out$ch)){
      warning(paste("\narguments warning: got dark channels '", unique(dark_channels$ch), 
                    "', but expected '", unique(arguments.df.out$ch), "'.",
                    "This may cause unexpected behaviour!\n",
                    collapse = " "))
    }
    # Add columnt to arguments
    arguments.df.out <- 
      dplyr::left_join(arguments.df.out,
                       dark_channels,
                       by="ch")
  }
  
  # List flat files
  flat_files <- dir(path, pattern = flat_pattern)
  if(length(flat_files) > 0){
    # Prepare DF
    flat_channels <- data.frame(flat = flat_files,
                                ch = sub(flat_pattern, "\\1", flat_files))
    # Compare correction channels with image channels 
    if(!setequal(flat_channels$ch, arguments.df.out$ch)){
      warning(paste("\narguments warning: got flat channels for'", unique(flat_channels$ch), 
                    "', but expected '", unique(arguments.df.out$ch), "'.",
                    "This may cause unexpected behaviour!\n",
                    collapse = " "))
    }
    # Add columnt to arguments
    arguments.df.out <- 
      dplyr::left_join(arguments.df.out,
                       flat_channels,
                       by="ch")
  }
  
  # Final checks ####
  if(all(is.na(arguments.df.out$t.frame))){
    warning("arguments warning: No t.frame data extracted, replacing all NAs with '0'. Check your directories and file.pattern if this is unexpected.")
    arguments.df.out$t.frame <- 0
  }
  
  if(any(is.na(arguments.df.out)) | any(arguments.df.out == "")){
    print(arguments.df.out)
    stop("arguments error: at least one of the values in the arguments.df dataframe is missing or blank, check your directories and file.pattern")
  }
  
  # Add Cell-ID's t.frame (annoying) indexing
  cids.frame <- arguments.df.out %>% 
    group_by(pos, t.frame) %>% 
    summarise(.groups = "drop_last") %>%
    arrange(pos, t.frame) %>% 
    mutate(t.frame.cid = 1:n() - 1)
  # And join it to the output df
  arguments <- left_join(arguments.df.out,
                         cids.frame, 
                         by = c("pos", "t.frame"))
  
  # Check for existing output files and stuff
  arguments_check(arguments, check_fail)
  
  # Chin pum!
  return(arguments)
}

#' Check if output files have already been created
#' 
#' @details 
#' 
#' Emits warnings, and optionally an error, to prevent accidental overwriting of existing output files.
#' 
#' Also checks if Cell-IDs t.frame numbering will match the t.frame in the file's names.
#' 
#' Also checks that the channel identifier is 3 characters long, as required by Cell-ID to distinguish imaging channels.
#'
#' @param arguments_df The table generated by \code{arguments}.
#' @param check_fail Whether to produce an error if any output file or directory already exists.
#'
#' @return TRUE if no existing output files were found.
#' @export
#'
arguments_check <- function(arguments_df, check_fail=F){
  n_row <- nrow(arguments_df)
  result <- TRUE
  
  # Check primary key ####
  # Check that t.frame indices from files will match Cell-ID's indices.
  t.frame.idx.match <- arguments_df %>% 
    select(t.frame, t.frame.cid) %>% unique() %>% 
    mutate(matches = t.frame == t.frame.cid)
  if(sum(!t.frame.idx.match$matches) > 0){
    warning(paste0(
      "\n\narguments_check: ",
      "Time frame index numbers in file names will not match those assigned by Cell-ID ",
      "in the 'out_all' files. This is likely to cause some confusion when filtering by 't.frame' ",
      "or when joining metadata to 'cdata'. Cell-ID's time frames start at 0 and grow by 1 unit.",
      "\n\n"
    ))
  }
  # Check channels ####
  # Check if the channel ID is of the right length
  ch_ids <- unique(arguments_df[["ch"]])
  ch_len_test <- any(nchar(ch_ids) != 3)
  if(ch_len_test){
    # Emit a warning
    warning(paste0(
      "\n\narguments warning: the following channel identifiers are of the incorrect length, which must be equal to 3 characters: '",
      paste0(ch_ids[nchar(ch_ids) != 3], collapse = "', '"),
      "'. Cell-ID will group imaging channels by the first 3 letters of file names. This is not customizable yet.\n\n"
    ))
    # And check if this will cause problems with Cell-ID, and stop if it is the case:
    unique_ch_check <- arguments_df %>% 
      mutate(first_ch_chars = substr(image, 1, 3)) %>% 
      select(first_ch_chars, ch) %>% 
      unique() %>% nrow()
    if(unique_ch_check != length(ch_ids)) {
      stop(paste(
        "\n\narguments error: the first three letters of fluorescence",
        "images and the extracted channel IDs form different sets.",
        "Rename your files and rebuild the arguments dataframe.\n\n"
      ))
    } else {
      warning(paste(
        "\n\narguments warning: though the channel IDs are short, there will be",
        "no issues with the current image set during imaging channel mapping.\n\n"
      ))
    }
      
  }
  
  # Check output files ####
  n_out <- sum(file.exists(arguments_df$output))
  if(n_out>0){
    warning(paste0("\n\narguments_check: ", 
                   n_out, "/", n_row, 
                   " output directories already exist.\n\n"))
  } else {
    message("\n\narguments_check: no output directories exist.\n\n")
  }
  
  # Check BF exists ####
  n_bf <- sum(file.exists(paste0(arguments_df$path, .Platform$file.sep, arguments_df$bf, ".out.tif")))
  if(n_bf>0){
    # Warning
    warning(paste0("\n\narguments_check: ", 
                   n_bf, "/", n_row, 
                   " output BF.out.tif files already exist.\n\n"))
  } else {
    # Message OK
    message("\n\narguments_check: no output BF.out.tif files exist.\n\n")
  }
  
  # Check FL exists ####
  n_fl <- sum(file.exists(paste0(arguments_df$path, .Platform$file.sep, arguments_df$image, ".out.tif")))
  if(n_fl>0){ 
    # Warning
    warning(paste0("\n\narguments_check: ", 
                   n_fl, "/", n_row, 
                   " output FL.out.tif files already exist.\n\n"))
  } else{
    # Message OK
    message("\n\narguments_check:  no output FL.out.tif files exist.\n\n")
  }
  
  
  # Raise an error if requested
  if (any(c(n_out, n_bf, n_fl) > 0)) {
    result <- FALSE
    if(check_fail) stop("\n\narguments_check: output files found, raising error!\n\n")
  }
  
  # Count channels ####
  ch_count <- arguments_df |> group_by(ch) |>
    summarise(image_count = n()) |> with(image_count) |> unique()
  if(length(ch_count) != 1){
    msg = paste0("\n\narguments_check: the amount of images differ between channels! (", 
                 paste(ch_count, collapse = " "),
                 ")\n\n")
    if(check_fail) {stop(msg)} else {warning(msg)}
  }
  
  # Count positions ####
  pos_count <- arguments_df |> group_by(pos) |>
    summarise(image_count = n()) |> with(image_count) |> unique()
  if(length(pos_count) != 1) {
    msg = paste0("\n\narguments_check: the amount of images differ between positions! (",
                 paste(pos_count, collapse = " "),
                 ")\n\n")
    if(check_fail) {stop(msg)} else {warning(msg)}
  }
  
  # Count frames ####
  t.frame_count <- arguments_df |> group_by(t.frame) |>
    summarise(image_count = n()) |> with(image_count) |> unique()
  if(length(t.frame_count) != 1) {
    msg = paste0("\n\narguments_check: the amount of images differ between t.frames! (",
                 paste(t.frame_count, collapse = " "),
                 ")\n\n")
    if(check_fail) {stop(msg)} else {warning(msg)}
  }
  
  return(result)
}

#' Default parameters list for Cell-ID
#' 
#' Returns a list of key-value pairs, for the default Cell-ID parameters.
#' It's output will tipically be used by \code{parameters_write}.
#' 
#' @details 
#' 
#' Documentation for each parameter can be found at: https://github.com/darksideoftheshmoo/cellID-linux#parameters
#' 
#' Boolean values are for "flag" type parameters
#' which enable a feature when present (eg. "align_fl_to_bf"),
#' or, if absent, indicate default behavior.
#' 
#' Other parameters have values
#' which must end up separated from names by a space " "
#' in the parameters.txt file format that Cell-ID uses:
#' 
#' \preformatted{
#' max_split_over_minor 0.5
#' max_dist_over_waist 8
#' max_pixels_per_cell 2000
#' min_pixels_per_cell 75
#' background_reject_factor 0.75
#' tracking_comparison 0.2
#' align_fl_to_bf
#' image_type brightfield
#' bf_fl_mapping list
#' }
#' 
#' @param max_split_over_minor Default: \code{0.50} For every combination of two pixels on the boundary, Cell-ID calculates the distance along the boundary path divided by the Euclidean distance between them. The maximum value of this ratio is larger for cells with a “figure-eight” shape that were pinched in some part than for circular cells. If the maximum value is above a user-defined threshold (which defaults to max_dist_over_waist=6), then the cell is split into two cells at the location of the pinch. After a split, if the Euclidean distance divided by the length of the minor axis of either of the new cells is greater than a user-defined value (which defaults to max_split_over_minor=0.5), then the two cells are re-grouped as a single cell. Thus, to perform the split we require that the two new cells have a generally circular shape and are not too elongated, as would be the case if the previous split was not over two cells, but over a cell and its mating projection.
#' @param max_dist_over_waist Default: \code{8.00} For every combination of two pixels on the boundary, Cell-ID calculates the distance along the boundary path divided by the Euclidean distance between them. The maximum value of this ratio is larger for cells with a “figure-eight” shape that were pinched in some part than for circular cells. If the maximum value is above a user-defined threshold (which defaults to max_dist_over_waist=6), then the cell is split into two cells at the location of the pinch. After a split, if the Euclidean distance divided by the length of the minor axis of either of the new cells is greater than a user-defined value (which defaults to max_split_over_minor=0.5), then the two cells are re-grouped as a single cell. Thus, to perform the split we require that the two new cells have a generally circular shape and are not too elongated, as would be the case if the previous split was not over two cells, but over a cell and its mating projection.
#' @param max_pixels_per_cell Default: \code{2000} Area limits per cell (upper bound, in pixels).
#' @param min_pixels_per_cell Default: \code{75} Area limits per cell (lower bound, in pixels).
#' @param background_reject_factor Default: \code{0.75} CellID's code makes an initial decision about the graylevels of the boundary pixels. To do this it takes the mean position of all the graylevels in the images and subtracts Z standard deviations. It then starts by considering all gray levels below this value as being parts of the cell borders. This value Z is the parameter background_reject_factor. Brightfield images taken slightly out of focus may do better with with higher values (ie, higher values will better avoid spurious cells), but if the cell boundaries in the image are too narrow, a smaller value may be necessary--which might increase the level of background.
#' @param tracking_comparison Default: \code{0.20} Cell-ID attempts to track cells over time. The value of this parameter is the minimal fractional overlap between two cells in consecutive time points for them to be considered the same cell. The default value is 0.2. Also named \"I_over_U_for_match\" in CellID's cell.c and segment.c files.
#' @param align_individual_cells Default: \code{F} Allow wiggling between the brightfield and fluorescence images.
#' @param align_fl_to_bf Default: \code{F} Frame alignment. Cell-ID can perform image registrations, moving the the frame in XY to align it to a reference image. If “align FL to BF” is selected the bright field image is used as reference. If “align FL to first” is selected the first fluorescence image is used as reference. These options are especially useful when sampling different positions, as repositioning of the microscope stage might introduce some displacement between consecutive images.
#' @param align_fl_to_first Default: \code{F} To-do: document or link to explanation.
#' @param image_type Default: \code{"brightfield"} To-do: document or link to explanation.
#' @param bf_fl_mapping Default: \code{"list"} Possible values: "list", "time". "bf_fl_mapping" option description (guessed from code, mask_mod branch). The mapping between brightfield and fluorescence images can be made by acquisition time, or derived from the order in the list of paths passed as command line options "-b" and "-f" to cell. If the order is by "list", then the paths must be grouped and ordered first by t.frame (ascending) and then by channel. If the order is by "time", cell derives the BF-FL mapping from the acquisition time in the TIFF metadata.
#' @param treat_brightfield_as_fluorescence_also Default: \code{F} Calculate all the fluorescence images variables on the bright field image as if it were a fluorescence image. This is potentially a good idea since it allows a good way to reject spurious cells. For example, the average value of the boundary pixels in good cells will be lower than the background level, but not so for spurious cells, etc.
#' @param third_image  Set to "vacuole_label" or "nuclear_label" to enable third image processing in Cell-ID. These images show the location of the nucleus or vacuole, and are used to derive measurements specific to those structures.
#' @return A nice list of named parameters, input for \code{parameters_write}.
#' 
#' @export
#' 
#' @seealso \link[rcell2.cellid]{parameters_write}, \link[rcell2.cellid]{arguments}
#' 
parameters_default <- function(
  max_split_over_minor = 0.50,
  max_dist_over_waist = 8.00,
  max_pixels_per_cell = 2000,
  min_pixels_per_cell = 75,
  background_reject_factor = 0.75,
  tracking_comparison = 0.20,
  align_individual_cells = F,
  align_fl_to_bf = F,
  align_fl_to_first = F,
  image_type = "brightfield",
  bf_fl_mapping = "list",
  treat_brightfield_as_fluorescence_also = F,
  third_image=F){
  
  return(list(
    max_split_over_minor = max_split_over_minor,
    max_dist_over_waist = max_dist_over_waist,
    max_pixels_per_cell = max_pixels_per_cell,
    min_pixels_per_cell = min_pixels_per_cell,
    background_reject_factor = background_reject_factor,
    tracking_comparison = tracking_comparison,
    align_individual_cells = align_individual_cells,
    align_fl_to_bf = align_fl_to_bf,
    align_fl_to_first = align_fl_to_first,
    image_type = image_type,
    bf_fl_mapping = bf_fl_mapping,
    treat_brightfield_as_fluorescence_also = treat_brightfield_as_fluorescence_also,
    third_image = third_image
  ))
}

#' Write parameters to a [temporary] file
#' 
#' Parses a \code{parameters.list} object from \code{parameters_default},
#' and writes its contents to a Cell-ID friendly plain text file.
#' 
#' @param parameters.list a parameters list for Cell-ID (like one from parameters_default)
#' @param param.dir directory where parameter files will be written.
#' @param param.file a file name for the parameters file.
#' @return A path to the text file where parameters where written.
#' 
#' @export
#' 
#' @seealso \link[rcell2.cellid]{parameters_default}, \link[rcell2.cellid]{arguments}
#' 
parameters_write <- function(parameters.list = rcell2.cellid::parameters_default(), 
                             param.dir = NULL,
                             param.file = NULL){
  
  # Print target directory if NULL
  if(is.null(param.dir)) param.dir <- base::tempdir()
  if(is.null(param.dir) & !is.null(param.file)) cat(paste0("\nSaving parameters file to: ", param.dir, "\n"))
  
  # Check if directory exists
  param.dir <- normalizePath(param.dir, mustWork = T)
  
  if(is.null(param.file)){
    param.file <- tempfile(tmpdir = param.dir, pattern = "parameters_", fileext = ".txt")
    cat(paste0("\nSaving parameters file to: ", param.file, "\n"))
  } else {
    param.file <- file.path(param.dir, param.file)
  }
  
  # Process the list into a valid parameter list
  # converting values to character type
  param_array <- 
    sapply(seq_along(parameters.list), function(i) {
      # For each parameter, get its name and value
      item_val <- parameters.list[[i]]
      item_name <- names(parameters.list)[i]
      
      # Boolean values are for "flag" type parameters
      # Which enable a feature when present (eg. "align_fl_to_bf")
      if(isTRUE(item_val))
        r <- item_name
      # And, if absent, indicate default behavior:
      else if(isFALSE(item_val))
        r <- NA
      # Other parameters have values
      # which must be separated from names by a space " "
      else
        r <- paste0(item_name, " ", item_val)
      
      # return "r" to sapply
      return(r)
    })
  
  # Filter any NAs (which come from FALSE flag-type parameters)
  param_array <- param_array[!is.na(param_array)]
  
  # Write to the parameter file
  write(x = param_array, file = param.file)
  
  return(param.file)
}


#' Load Cell-ID's output files into R
#'
#' @param path Path to Cell-ID's output directory, tipically also the images directory.
#' @param pdata Path to metadata CSV file.
#' @param position.pattern Regex describing what the position string looks like (default ".*Position(\\d+).*") including a capturing group for the position ID number (coerced to integer).
#' @param fluorescence.pattern Regex describing what the fluorescence/channel ID string looks like (default "^([GCYRT]FP|[GCYRT]\\d+)_Position\\d+_time\\d+.tif$"). There must be only one capturing group, ant it must be for the channel identifier.
#' @param ucid.zero.pad Amount of decimal digits for the cellID (defaults 4, corresponding to a maximum of 9.999 cellIDs and 9999 positions).
#' @param append.posfix String appended to the channel ID extracted by `fluorescence.pattern` (`NULL` by default, but "FP" is usual).
#' @param ... Arguments passed on to \code{load_out_all}. Patterns for "out" files, fluorescence channel, and other options may be changed here.
#' @inheritDotParams load_out_all
#' @return A list of dataframes: data (CellID data), images (images metadata and paths), image_maping (extra mapping metadata from CellID: BF to FL correspondence, channel flag, bf_as_fl flag, and one-letter channel encoding).
# @examples
# cell_data <- get_cell_data(path = path, pdata = pdata)
#' @import dplyr stringr tidyr readr
#' @importFrom purrr map
#' @export
get_cell_data <- function(path,
                          pdata = NULL,
                          position.pattern = ".*Position(\\d+).*",
                          # fluorescence.pattern = "^([GCYRT]FP)_Position\\d+.tif$",
                          # fluorescence.pattern = "^(BF|[GCYRT]FP|[GCYRT]\\d+)_Position\\d+_time\\d+.tif$",
                          fluorescence.pattern = "^(BF|[GCYRT]FP|[GCYRT]\\d+)_Position\\d+.*.tif$",
                          ucid.zero.pad = 4,
                          append.posfix = NULL,
                          ...){
  # Normalize the path
  path <- normalizePath(path, mustWork = T)
  
  # Load data ####
  # Cargar datos out_all y juntar en un solo dataframe usando metadata de "out_bf_fl_mapping"
  cat("\n\nLoading CellID output files...\n")
  d.list <- load_out_all(path = path,
                         position.pattern = position.pattern,
                         fluorescence.pattern = fluorescence.pattern,
                         ...)  # https://stackoverflow.com/questions/40794780/r-functions-passing-arguments-with-ellipsis/40794874
  
  # Create ucid column ####
  cat("\rCreating ucid column...                            ")
  d.list$d <- d.list$d %>%
    # By padding to an invariant number of digits, cells from positions as 90 and 9 _could_ have the same UCID.
    # The next lines fix that bug, which would cells to not be filtered correctly, and be plotted anyways, or other problems.
    # Also, the padding should not be very large, or it will "overflow" R's integer class.
    # To-do: replace the following code with str_pad
    mutate(cellid.pad = ucid.zero.pad - nchar(as.character(cellID))) %>%
    mutate(
      ucid = as.integer(
        paste0(pos,
               sapply(cellid.pad, FUN = function(pad) paste0(rep(x="0", 
                                                                 times=max(0, pad)), 
                                                             collapse = "")),
               cellID))) %>% 
    select(ucid, tidyselect::everything())
  
  if(any(d.list$d$cellid.pad < 0))
    stop("cellID too large to pad, increase ucid.zero.pad (and check for integer overflow).")

  # Delete the cellid.pad column
  d.list$d$cellid.pad <- NULL

  # Merge with pdata ####
  # if(exists("pdata", inherits = F)){
  cat("\rJoining pdata if specified...")
  if(!is.null(pdata)){
    pdata <- readr::read_csv(pdata) %>% mutate(pos = pos %>% as.numeric)
    d.list$d <- d.list$d %>% left_join(pdata, by = "pos")
    cat(" and it was :)                            ")
  } else cat(" but it was not :(                            ")

  # Create paths DF ####
  # Create paths dataframe and add three-letter code for channel
  cat("\rCreating image paths dataframe...                            ")
  paths <- d.list$d.map
  if(!is.null(append.posfix)){
    paths <- mutate(paths, channel = paste0(toupper(channel), append.posfix))
  }

  # Bind paths dataframe it with itself, to get entries for BF as well.
  paths <- dplyr::bind_rows(
    paths %>%  # Get FL image paths
      select(pos, t.frame, channel, fluor) %>%
      rename(file = fluor),
    
    paths %>%  # Get BF image paths
      select(pos, t.frame, channel, bright) %>%
      rename(file = bright) %>%
      mutate(channel = "BF") %>% 
      unique()  # There may be BF path duplicates in the "bright" column, so keep the unique set
  ) %>%
    mutate(path = dirname(file),  # Add the directory path
           is.out = FALSE)        # Add the is.out column

  paths <- bind_rows(paths,
                     paths %>%  # bind it with part of itself, out files are named exactly the same, but with an extra ".out.tif"
                       mutate(file = paste0(file, ".out.tif"),
                              channel = paste0(channel, ".out"),
                              is.out = TRUE)) %>%
    mutate(image = basename(file))
  # Save to the output list
  d.list$d.paths <- paths
  
  # Position directory df ####
  positions <- data.frame(
    pos = as.numeric(str_replace(string = basename(d.list$pos.directories), pattern = position.pattern, replacement = "\\1")),
    output = d.list$pos.directories
  )
  
  # Check that image paths exists, or try to fix them using the path to the data.
  images <- check_and_fix_paths(path=path, images=d.list$d.paths)
  
  # Make output list ####
  cell.data <- list(data = d.list$d,
                    images = images,
                    mapping = d.list$d.map,
                    channels = unique(d.list$flag.channel.mapping),
                    positions = positions,
                    variable_descriptions = cellid_output_descriptions())
  
  cat("\rDone loading CellID data!                            \n")
  
  # Return ####
  return(cell.data)
}

#' @rdname get_cell_data
#' @export
cell.load.alt <- get_cell_data

#' Check and fix image paths
#' 
#' @param path Path to the data directory, holding the images.
#' @param images Cell-ID dataframe with image paths, as loaded by \code{get_cell_data}.
#' 
#' @export
check_and_fix_paths <- function(path, images){
  # TEST:
  # images <- d.list$d.paths
  
  if(any(!file.exists(images$file))){
    warning("\nNot all image files exist in the expected filesystem directory. Attempting to fix them... ")
    new_paths <- file.path(path, basename(images$file))
    if(all(file.exists(new_paths))){
      warning("Image paths fixed.\n")
      images$path <- path
      images$file <- file.path(path, basename(images$file))
    } else {
      warning("Could not fix image file paths, expect issues when loading images in rcell2.\n")
    }
  }
  
  return(images)
}

#' Una función que lea un .csv y les agregue una columna con un id del archivo (pos)
#' @keywords internal
#' @import stringr dplyr
read_tsv.con.pos <- function(.path.archivo, 
                             position.pattern, 
                             col_types = "c"){
  # Intro message:
  cat(paste0("\rReading: '", basename(.path.archivo), 
             "' in directory '", basename(dirname(.path.archivo)), "'.") ,
      "\033[K")
  
  # Test:
  # .path.archivo <- .nombre.archivos[1]
  
  # Normalize file name
  .archivo <- normalizePath(.path.archivo, mustWork = T)
  # Get position index number
  .pos <- stringr::str_replace(dirname(.path.archivo), position.pattern, "\\1") %>% 
    as.integer()
  
  # Load "out_all" files ####
  d <-  readr::read_tsv(.archivo, col_types = col_types, trim_ws = T) %>%
    # "pos", "t.frame" y "flag" estan bf_fl_mapping y en out_all
    mutate(pos = as.integer(.pos),  # La columna de ID es "pos"
           t.frame = as.integer(t.frame),
           flag = as.integer(flag))
    # el resto de las columnas que no se comparten y deberian ser enteras
    # se convierten en load_out_all()
  
  # Check "con.vol" duplicate column name (old CellID bug)
  if("con.vol_1" %in% names(d)) {
    cat(paste0("\nRemoving 'con.vol_1' column from position: ", 
               .pos, 
               ". Use CellID version > 1.4.6 to stop seeing this message.\n"))
    d <- dplyr::select(d, -con.vol_1)
  }
  con.vol.dupes <- startsWith(names(d), "con.vol")
  if(sum(con.vol.dupes) > 1){
    first <- names(d)[which(con.vol.dupes)][1]
    dupes <- names(d)[which(con.vol.dupes)][-1]
    dupes <- paste(dupes, collapse = ", ")
    
    cat(paste0("\nRemoving '", dupes, "' column(s) from position: ", 
               .pos, 
               ". Use CellID version > 1.4.6 to stop seeing this message.\n"))
    
    colnames(d)[which(con.vol.dupes)[1]] <- "con.vol"
    d[which(con.vol.dupes)[-1]] <- NULL
  }
  
  # Done
  return(d)
}

#' Una función para leer y combinar todos los archivos "out".
#'
#' @inheritParams get_cell_data
#' @param out_file_pattern Regex matching CellID's main output file.
#' @param out_mapping_pattern Regex matching CellID's image mapping output file.
#' @param id_columns_check Whether to check for problems in "ID" and "value" columns. If TRUE, the function checks the uniqueness of ID columns.
#' @param fix_id_columns Whether apply the \code{fix_id_columns_fun} summary function to problematic ID columns (with non-unique vales per cell across channels).
#' @param fix_id_columns_fun Summary function for \code{fix_id_columns}, uses dplyr's \code{first} by default (thereby using the value of the first channel).
#' @import dplyr tidyr readr
#' @importFrom purrr map
#' @importFrom data.table setDT dcast setDF
#' @return A list of two dataframes: `d` contains the actual output, and `out.map` contains image paths and metadata.
#' @keywords internal
load_out_all <- function(path,
                         position.pattern = ".*Position(\\d+)$",
                         out_file_pattern = "^out_all$",
                         out_mapping_pattern = "^out_bf_fl_mapping$",
                         fluorescence.pattern = ".*(BF|[A-Z]FP)_Position.*",
                         id_columns_check=TRUE,
                         fix_id_columns=FALSE,
                         fix_id_columns_fun=dplyr::first){
  
  # List position directories
  .pos.directories <- list.dirs(path = path, full.names = T) %>% 
    grep(pattern = position.pattern, value = T)
  
  # List paths de los "out_all"
  .nombre.archivos <- list.files(path = .pos.directories, full.names = T,
                                 pattern = out_file_pattern)
  
  # List paths de los "bf_fl_mapping"
  .nombre.archivos.map <- list.files(path = .pos.directories, full.names = T,
                                     pattern = out_mapping_pattern)
  # A bit of error handling
  if(length(.nombre.archivos) == 0) 
    stop("Error in load_out_all: no CellID output files found, check your path, options and files.")
  if(length(.nombre.archivos.map) == 0) 
    stop("Error in load_out_all: no CellID mapping files found, check your path, options and files.")
  if(length(.nombre.archivos) != length(.nombre.archivos.map)) 
    stop("Error in load_out_all: different amount of mapping and cellid output files.")
  
  # Cargo y junto los "out_all"
  cat("\rLoading datasets...\033[K")
  d.out <- purrr::map(.x = .nombre.archivos,
                      .f = read_tsv.con.pos, 
                      # .carpeta = path,  # Argumento reemplazafo por ".path.archivo".
                      # col_types = "iiiiiddddddddddddddddddddddddddddddddddddddddddddddddddd", # types: 5 int columns, 51 double columns
                      col_types = readr::cols(.default = "d"), # types: all double, convert later
                      position.pattern = position.pattern) %>%
    bind_rows() %>% 
    mutate(cellID = as.integer(cellID))
  
  cat("\n Done loading 'out_all' files!\n")

  # # Cargo y junto los "out_bf_fl_mapping"
  cat("\rLoading mapping...              ")
  d.map <- purrr::map(.x = .nombre.archivos.map,
                      .f = read_tsv.con.pos,  # Una función para leer los archivos "out" y agregarles "pos" segun la carpeta que los contiene
                      # .carpeta = path,  # Argumento reemplazafo por ".path.archivo".
                      col_types = "ciicl",  # types: 2 char columns, 2 int columns, 1 logical column
                      position.pattern = position.pattern) %>%
    bind_rows() %>%
    mutate(channel = str_replace(string = basename(fluor),
                                 pattern = fluorescence.pattern,
                                 replacement = "\\1")) %>%
    mutate(channel = tolower(channel))
  
  # Check if any fluorescence images are missing
  missing_fl_test <- d.map %>% select(pos, flag, t.frame) %>% unique() %>% 
    group_by(pos, t.frame) %>% summarise(count = n(), .groups = "drop") %>% 
    with(length(unique(count)) == 1)
  if(!missing_fl_test){
    warning("One or more fluorescence images are missing from the set. This will generate missing values or cause unexpected errors.\n")
  }
  
  # keep flag-channel mapping for later...
  flag.channel.mapping <- unique(dplyr::select(d.map, flag, channel))
  
  # Check amount of flags VS amount of extracted channel names
  flag.ch.check <- length(unique(flag.channel.mapping$flag)) == length(unique(flag.channel.mapping$channel))
  if(!flag.ch.check) 
    stop(
      "Error in load_out_all: number of flags does not match extracted channel names.\nCheck your 'fluorescence.pattern'\n\n",
      paste0("Flags: \n", paste(unique(flag.channel.mapping$flag), collapse = "\n"), "\n\n"),
      paste0("Channels: \n", paste(unique(flag.channel.mapping$channel), collapse = "\n"), "\n\n")
    )
  
  cat("\n Done loading 'bf_fl_mapping' files!\n")
  avail.ch <- paste(flag.channel.mapping$channel, flag.channel.mapping$flag, sep = "=", collapse = ", ")
  cat(paste("    Final channel-flag mapping:", avail.ch, "\n"))

  # Return join (discard flag variable)
  cat("\rJoining data and mapping...\033[K")
  d.out.map <- dplyr::left_join(
      # Data:
      d.out,
      # Image mapping:
      unique(select(d.map, flag, t.frame , pos, channel)),
      by = c("flag", "t.frame", "pos")
    ) %>%
    select(-flag)
  
  # Check row numbers
  if(nrow(d.out.map) > nrow(d.out)){
    stop("Error at load_out_all: while joining output and mapping, at least one output row matched multiple mappings.")
  }
  
  # Data table leftjoin tests
  if(F){
    mapping <- unique(select(d.map, flag, t.frame , pos, channel))
    data.table::setDT(mapping)
    data.table::setDT(d.out)
    r <- d.out[mapping, on = c("flag", "t.frame", "pos")]
  }
  
  # Add f.tot columns to data
  cat("\rCreating FL variables...                           ")
  d.out.map <- mutate(d.out.map,
                      f = f.tot - (a.tot * f.bg),
                      cf = f / a.tot,
                      f.loc = f.tot - (f.local.bg * a.tot),
                      cf.loc = f.loc / a.tot)
  
  
  # Calculate el.p ####
  ellipse_perimeter <- function(maj.axis, min.axis){
    pi * (3 * (maj.axis / 2 + min.axis / 2) - sqrt((3 * maj.axis / 2 + min.axis / 2) * (maj.axis / 2 + 3 * min.axis / 2)))
  }
  # ellipse.perim = perimeter of theoretical ellipse, calculated using each
  # cell's axis values.
  # el.p = ratio of ellipse perim over the perimeter measured by cellID.
  # If this number is small ( < ~0.7) it's probably not a cell.
  cat("\rCreating el.p column...                            ")
  d.out.map <- dplyr::mutate(d.out.map,
                             ellipse.perim = ellipse_perimeter(maj.axis, min.axis),
                             el.p = ellipse.perim / perim)
  
  # d.out.map <- filter(d.out.map, cellID==1)  # test one cell
  
  # Define ID and value variables ####
  # Variables that should not change across fluroescence channels are used as IDs
  cell_idcols <- c("cellID", "t.frame", "pos")
  id_cols = c("cellID",
              # flag,  # dropped above
              "t.frame",
              "time",
              "pos",
              "xpos",
              "ypos",
              "a.tot",
              "num.pix",
              "fft.stat",
              "perim",
              "maj.axis",
              "min.axis",
              "ellipse.perim",
              "el.p",
              "rot.vol",
              "con.vol",
              "a.tot.p1",  # should be moved from the id_cols, issue #32
              "a.tot.m1",
              "a.tot.m2",
              "a.tot.m3",
              # a.local.bg,  # moved from the id_cols, issue #32
              "a.local",
              # a.local2.bg,  # moved from the id_cols, issue #32
              "a.local2",
              "a.surf",
              # con.vol_1,  # duplicated, removed by read_tsv.con.pos and in recent CellID versions
              "sphere.vol")

  id_cols_notcell <- setdiff(id_cols, cell_idcols)
  
  values_from = c("f.tot",
                  "f", 
                  "cf",
                  "f.loc",
                  "cf.loc",
                  "f.nucl",
                  "a.nucl",
                  "a.vacuole",
                  "f.vacuole",
                  "a.local.bg", # moved from the id_cols, issue #32
                  "a.local2.bg",  # moved from the id_cols, issue #32
                  "f.bg",
                  "f.tot.p1",
                  "f.tot.m1",
                  "f.tot.m2",
                  "f.tot.m3",
                  "xpos.nucl",
                  "ypos.nucl",
                  "f.nucl1",
                  "f.nucl.tag1",
                  "a.nucl1",
                  "f.nucl2",
                  "f.nucl.tag2",
                  "a.nucl2",
                  "f.nucl3",
                  "f.nucl.tag3",
                  "a.nucl3",
                  "f.nucl4",
                  "f.nucl.tag4",
                  "a.nucl4",
                  "f.nucl5",
                  "f.nucl.tag5",
                  "a.nucl5",
                  "f.nucl6",
                  "f.nucl.tag6",
                  "a.nucl6",
                  "f.local.bg",
                  "f.local2.bg")
  
  cat("\rChecking ID column uniqueness...\033[K")

  # ID cols checks ####
  # Check if "id columns" are really the same within each observation
  if(id_columns_check){
    id_cols_check <- d.out.map[,id_cols] %>% group_by_at(.vars = cell_idcols) %>% 
      summarise_all(.funs = function(x) length(unique(x)) == 1) %>% 
      arrange(pos, cellID, t.frame)
    
    id_cols_check$any_bad <- !apply(id_cols_check[,id_cols_notcell], 
                                    MARGIN = 1, 
                                    FUN = all)  
  
    if(any(id_cols_check$any_bad)){
      
      which_bad <- id_cols_check %>%
        filter(any_bad) %>% .[,id_cols_notcell] %>% {split(., 1:nrow(.))} %>% 
        lapply(FUN = function(x) {
          which_bad_idx <- which(!x)
          
          # Handle all-good columns with NA value (shouldnt be necessary but whatever).
          if (length(which_bad_idx) == 0) {
            which_bad <- c(NA_character_)
          } else {
            which_bad <- c(id_cols_notcell[which_bad_idx])
          }
          
          return(which_bad)
        })
      
      # Unlist!
      which_bad <- which_bad %>% unlist() %>% unique()
      
      # Fix the problem
      if (length(which_bad) > 0) {
        # Prepare warning text
        warning_text <- paste(
          "Columns [", paste(which_bad, collapse = ", "), "] are not ID columns,",
          "they have different values across channels of the same cell.",
          "This might be due to enabling 'align_fl_to_bf' in Cell-ID."
        )
        
        # Choose how to handle the problem
        if(!fix_id_columns){
          # Warn the user
          warning(paste(warning_text, 
                        "Treating these variables as 'value' columns.",
                        "The column names in cdata will change!",
                        "Try enabling 'fix_id_columns' and setting an appropriate 'fix_id_columns_fun' argument.\n"))
          # Remove "bad" columns from the "id" vector
          id_cols <- id_cols[!id_cols %in% which_bad]
          # And add "bad" columns to the "values" vector
          values_from <- unique(c(values_from, which_bad))
        } else {
          # Warn the user
          warning(paste(warning_text, 
                        "Applying the 'fix_id_columns_fun' summary function to the problematic columns.",
                        "Try disabling 'fix_id_columns' to treat these columns as regular 'value' columns.\n"))
          # Use the "fix_id_columns_fun" summary function to "fix" the ID variables with non-unique values
          d.out.map <- d.out.map %>% group_by_at(.vars = cell_idcols) %>% 
            arrange(channel) %>% 
            mutate_at(.vars = id_cols, .funs = fix_id_columns_fun) %>% ungroup()
        }
      }
    }
  }
  # Convert to wide format ####
  # Right now the out_all is in a "long" format for the "channel" variable.
  # Spread it to match expectations:
  cat("\rSpreading data from channels...\033[K")

  # Replacing pivot_wider with dcast (it is 2x faster)
  # cdata <- d.out.map %>%
  #   tidyr::pivot_wider(
  #     # Names for fluorescence columns come from the "channel" variable
  #     names_from = channel, names_sep = ".",
  #     id_cols = all_of(id_cols),
  #     values_from = all_of(values_from)
  #     )
  
  # New long to wide reshaping, with data.table:
  dcast_formula_lhs <- paste(collapse = " + ", id_cols)  # Using "cell_idcols" instead of "id_cols" is much faster, but the omitted columns would be lost
  dcast_formula_rhs <- "channel"
  dcast_formula <- paste0(dcast_formula_lhs, " ~ ", dcast_formula_rhs)
  dcast_values <- values_from
  # Reshape:
  data.table::setDT(d.out.map)
  cdata <- data.table::dcast(
    data = d.out.map, 
    formula = as.formula(dcast_formula), 
    value.var = dcast_values, sep = "."
  )
  data.table::setDF(cdata)
  
  # Check row number against expectation
  check_row_n <- nrow(cdata) == nrow(d.out.map)/length(unique(d.out.map$channel))
  if(!check_row_n) warning("Casting fluorescence columns produced the wrong number of rows per fluorescence channel.")
  
  cat("\rPreparing output...\033[K")
  # Unnecesary with DT dcast
  # cdata <- mutate(cdata, 
  #                 cellID = as.integer(cellID),
  #                 t.frame = as.integer(t.frame))

  # Prepare output list ####
  d.list <- list(
    "d" = cdata,
    "d.map" = d.map,
    "flag.channel.mapping" = flag.channel.mapping,
    "pos.directories" = .pos.directories
    )
  
  # Check unique ucid-frame combo ####
  # Check uniqueness of ucid-t.frame combinations
  if(nrow(unique(cdata[,c("cellID", "pos", "t.frame")])) < nrow(cdata)){
    
    dump.file <- tempfile(fileext = ".RDS")
    cat("\rFailed ucid-t.frame uniqueness check, dumping to RDS file:", dump.file ,"\033[K")
    saveRDS(d.list, dump.file)
    
    # Find problematic positions
    test.df <- cdata[,c("cellID", "pos", "t.frame")]
    test.df.list <- split(test.df, test.df$pos)
    test.res <- 
      lapply(test.df.list, function(d){
        nrow(unique(d)) < nrow(d)
      }) %>% unlist()
    
    stop(paste(
      "\nERROR: There are repeated cellID's in the out_all file! Dumped data to:",
      dump.file,
      "Problematic positions:", paste(names(test.res[test.res]), collapse = " ")
      )
    )
  }
  
  # Return ####
  return(d.list)
}

#' Print rcellid arguments summaries
#' 
#' A function to print some summaries, to check cellArgs2 output.
#' 
#' @param arguments The "arguments" dataframe, output from \code{rcell2.cellid::arguments()}.
#' @import dplyr
#' @export
arguments_summary <- function(arguments){
  arguments %>% group_by(ch) %>% summarise(n_count = n(), .groups = "drop") %>% print()
  arguments %>% select(bf) %>% summarise(unique_BF = "", n_count = length(unique(bf)), .groups = "drop") %>% print()
  arguments %>% group_by(t.frame) %>% summarise(n_count = n(), .groups = "drop") %>% print()
  arguments %>% group_by(pos) %>% summarise(n_count = n(), .groups = "drop") %>% print()
}

#' Make and "images" dataframe from "arguments" dataframe
#' 
#' The images dataframe is needed by many rcell2 functions. If it is not available from the output of \code{load_cell_data} or \code{get_cell_data}, then this function can help.
#' 
#' It essentially does a pivot_longer of the arguments.
#' 
#' @param arguments The "arguments" dataframe, output from \code{rcell2.cellid::arguments()}.
#' 
#' @return A data.frame similar to \code{get_cell_data()$images}.
#' @import dplyr
#' @export
arguments_to_images <- function(arguments){
  
  # BF dataframe
  bf <- arguments %>% 
    select(pos, t.frame, path, bf) %>% 
    unique() %>%
    mutate(ch = "BF") %>% 
    rename(image = bf)
  
  # FL dataframe
  fl <- arguments %>% select(pos, t.frame, path, image, ch)
  
  images <- 
    dplyr::bind_rows(bf, fl) %>% 
    mutate(file = file.path(path, image),
           is.out = F) %>% 
    select(pos, t.frame, ch, file, path, is.out, image) %>% 
    rename(channel = ch) #%>% 
  
  return(images)
}

#' Pipe
#'
#' purrr's pipe operator
#'
#' @importFrom purrr %>%
#' @name %>%
#' @rdname pipe
#' @param lhs,rhs specify what lhs and rhs are
# @examples
#' @keywords internal
#' # some examples if you want to highlight the usage in the package
NULL

#' Pad numbers to the same length with leading zeros
cero_a_la_izquierda <- function(x, pad_char="0"){
  set_pad <- ceiling(log10(max(x)+1))
  stringr::str_pad(string = x, width = set_pad, pad = pad_char)
}

#' Image file renamer for Metamorph MDA
#' 
#' MDA: "Multi dimensional acquisition" app in Metamorph.
#' 
#' Uses regex groups to extract channel, position and time information from file names, and uses it to stitch new and friendlyer names.
#' These are used to copy or link image files to a target directory.
#' 
#' For example, \code{far1_rtcc_exp16_thumb_w1LED-BF--YFPcube--cam_s17_t35.TIF} can be converted to \code{BF_Position17_time35.tif}.
#' 
#' The \code{identifier.pattern} is a key parameter. There must be three groups, one for each of the three information types: channel, position and time.
#' The defaults are useful for a file name such as "\code{far1_rtcc_exp16_thumb_w1LED-BF--YFPcube--cam_s17_t35.TIF}", in which the channel is identified by a "w", 
#' position by an "s", and time by a "t".
#' 
#' The order in which this information appears in the file name is specified in \code{identifier.info}. 
#' If you wish to add a prefix to each field in the final file name, name the elements in this vector. 
#' For example, the default \code{c("ch", Position="pos", time="t.frame")} indicates that channel has no prefix,
#' the "pos" field will be prefixed by "Position", and the "t.frame" field will be prefixed by "time".
#' Then, for example, a new file name could look like this: \code{BF_Position1_time3.tif}.
#' 
#' Channel names will be translated according to the rows in \code{channel.dict} (see the parameter's description).
#' These are easily adaptable to other use cases, for example you may change \code{channel.dict} to include more, less or other channels, in whatever order.
#' Note that the values in the \code{ch} column must exactly match the strings captured by the corresponding capture group in \code{identifier.pattern}.
#' For example, the channel in the original file names may be integers from 1 to 3, which are captured and matched with dplyr's left_join to the \code{channel.dict} data frame.
#' Then, the value in \code{ch.name} is used to build the final file name.
#' 
#' **Limitations**: In the original file names, the identifiers for each field can only be integers.
#'
#' examples 
#' images.path <- "~/Projects/PhD/data/uscope/multidimensional_exp-20211126-Far1NG-wt_y_dKar4/"
#' rename_mda(images.path, rename.function = file.copy)
#' 
#' @import dplyr
#' @param images.path Path to the directory containing the original images. Can be NULL if \code{file.names} is provided.
#' @param rename.path Path to the target directory. If \code{NULL} (the default) images are sent to a new "renamed" sub directory of \code{images.path}. The directory will be created if it does not exist.
#' @param rename.dataframe If \code{TRUE}, no renaming takes place and a data frame is returned instead. It can be altered by the user, and then passed to this same argument to rename the images with the adjusted names or mappings. The data frame is normally returned in the "status" item of the regular output.
#' @param rename.function Either \code{\link[base]{file.copy}}, \link[base]{file.symlink} or a similar function. Set to \code{NULL} to disable renaming (i.e. for testing purposes).
#' @param identifier.pattern Regex defining the iamge file pattern, with gropus for identifier in the file names.
#' @param identifier.info Character vector with strings "pos", "t.frame", and "ch" (channel), in the same order in which they appear in the \code{identifier.pattern}. If an element in the vector is named, the name is prefixed to the identifier in the final file name (for example, by default, "Position" is prepended to the position number; but channel has no prefix).
#' @param channel.maping.df A dataframe with two columns: "ch" holding the original channel names in the source files (coerced to character), and "ch.name" with the new names for each channel.
#' @param file.ext File extension to use in the final file name, such as: ".tif".
#' @param skip.thumbs.pat A regex pattern to filter out files. Convenient if the MDA output thumbnails for each image. Set to \code{NULL} to disable.
#' @param cleanup.first Set to \code{TRUE} to remove all files within the \code{rename.path} directory before renaming. \code{FALSE} by default.
#' @param file.names A character vector of image names to be processed; as an alternative to listing files in \code{images.path}. Ignored if \code{images.path} is provided.
#' @param ... Further arguments passed on to \code{rename.function}.
#' @export
#' @return Invisibly returns a list with the unique rename.paths (output directories), and a \code{data.frame} with the mappings for name conversions and the output value from the renaming function (see the \code{rename.function} parameter's description).
#' @import stringr dplyr
rename_mda <- function(images.path = NULL,
                       rename.path = NULL,
                       rename.dataframe = NULL,
                       rename.function = file.symlink,
                       identifier.pattern=".*_w(\\d).*_s(\\d{1,2})_t(\\d{1,2}).TIF$",
                       identifier.info = c("ch", Position="pos", time="t.frame"),
                       channel.maping.df = data.frame(ch=1:3, ch.name=c("BF", "YFP", "TFP")),
                       file.ext=".tif",
                       skip.thumbs.pat = ".*thumb.*",
                       cleanup.first=FALSE,
                       file.names = NULL,
                       ...
                       ){
  
  if(F){
    file.names <- c(
      "timecourse1_w1LED-BF--YFPcube--cam_s1_t1.TIF",
      "timecourse1_w1LED-BF--YFPcube--cam_s1_t10.TIF", 
      "timecourse1_w1LED-BF--YFPcube--cam_s1_t11.TIF")
    images.path = NULL
    rename.path = NULL
    rename.function = file.symlink
    identifier.pattern=".*_w(\\d).*_s(\\d{1,2})_t(\\d{1,2}).TIF$"
    identifier.info = c("ch", Position="pos", time="t.frame")
    channel.maping.df = data.frame(ch=1:3, ch.name=c("BF", "YFP", "TFP"))
    file.ext=".tif"
    skip.thumbs.pat = ".*thumb.*"
    cleanup.first=F
  }
  
  # Prepare the renaming data.frame unless one has been provided.
  if(!is.data.frame(rename.dataframe)){
    # Checks
    if(is.null(images.path) & is.null(file.names) & is.null(rename.dataframe)) {
      stop("rename_mda: error, either 'images.path', 'rename.dataframe' or 'file.names' must be specified.")}
    if(!( setequal(identifier.info, c("ch", "pos", "t.frame")) && length(identifier.info) == 3 )){
      stop("rename_mda: error, malformed identifier.info.")}
    
    # Get file names
    if(!is.null(images.path)){
      image.files <- dir(images.path, pattern = identifier.pattern, full.names = T)
    } else {
      image.files <- file.names
    }  
    
    # Skip thumbnails
    if(!is.null(skip.thumbs.pat)) image.files <- image.files[!grepl(skip.thumbs.pat, basename(image.files))]
    
    # Check if we got something
    if(length(image.files) == 0){
      stop(paste0(
        "rename_mda: error, no file names were retrieved from directory ",
        "'", images.path, "' ",
        "using pattern: '", gsub("([\\])","\\\\\\\\", identifier.pattern), "'"
      ))
    }
    
    # Extract groups
    images.info <- stringr::str_match(basename(image.files), identifier.pattern)[,-1,drop=F]
    # Convert to data frame
    images.info <- setNames(data.frame(images.info), identifier.info)
    
    # Checks
    channel.maping.df$ch <- as.character(channel.maping.df$ch)
    mapping_check <- all(unique(images.info$ch) %in% channel.maping.df$ch)
    if(!mapping_check) stop(paste0("rename_mda: error. Expected channel indexes '",
                                   paste(channel.maping.df$ch, collapse = ", "),
                                   "' and found channel indexes '",
                                   paste(unique(images.info$ch), collapse = ", "),
                                   "'. Consider updating 'channel.maping.df' to add a missing channel index."
    ))
    
    # Join extracted image info to the channel mapping data.frame
    images.info <- dplyr::left_join(images.info, channel.maping.df, by = "ch")
    
    # Add file name and path
    images.info$path <- images.path
    images.info$file <- image.files |> basename()
    
    # Figure out the path where renamed images will be written.
    if(is.null(rename.path)){ 
      if(is.null(images.path)){
        # If no "images.path" nor "rename.path" were specified, set "rename.path" to an empty string,
        # which ends up being the current working directory. Cleanup will be skipped in this case.
        rename.path <- ""
      } else{
        # If "rename.path" was not specified, use "renamed" as a sub-directory of "images.path".
        rename.path <- paste0(images.path, "/renamed")
      }
    }
    
    # Make new names and paths
    images.info$rename.path <- rename.path
    
  } else {
    warning("rename_mda: a data.frame was passed to 'rename.dataframe'.",
            " It's 'rename.file' column will be regenerated (overwritten).",
            " Will also update the 'rename.path' column if the parameter was specified.")
    
    # Update rename path if specified.
    images.info <- rename.dataframe
    if(!is.null(rename.path)){
      images.info$rename.path <- rename.path
    } else {
      rename.path <- images.info$rename.path |> unique()
    }
  }
  
  # Create the directory for renamed images.
  if(!is.null(rename.path)) dir.create(rename.path, showWarnings = F)
  
  # Delete images in "rename.path" if requested (unless it is empty).
  if(cleanup.first){
    if((rename.path != "") & (!is.null(rename.path))){
      # Delete all files recursively.
      cleanup.files <- dir(path = rename.path, full.names = T,
                           pattern = paste0(".*", file.ext))
      unlink(cleanup.files)
      
      # Re create output directory.
      dir.create(rename.path, showWarnings = F, recursive = T)
    } else {
      message("rename_mda: to prevent accidental deletions, cleanup was skipped because 'rename.path' is empty.")
    }
  }
  
  # Make new names.
  images.info$rename.file <- paste0(
    
    # ifelse(nzchar(names(identifier.info[identifier.info=="ch"])), "_", ""),
    names(identifier.info[identifier.info=="ch"]),
    images.info$ch.name,
    "_",
    
    names(identifier.info[identifier.info=="pos"]),
    cero_a_la_izquierda(as.numeric(images.info$pos)),
    "_",
    
    names(identifier.info[identifier.info=="t.frame"]),
    cero_a_la_izquierda(as.numeric(images.info$t.frame)),
    
    file.ext
  )
  
  # If a "renaming data frame" (images.info) was not provided, but "rename.dataframe" is TRUE,
  # then immediately return the prepared 'images.info' data frame.
  if(isTRUE(rename.dataframe)){
    message("rename_mda: rename.dataframe was set to TRUE, returning the 'images.info' dataframe. No renaming took place.")
    return(images.info)
  }
  
  # Rename the files.
  if(!is.null(rename.function)){
    status <- rename.function(from = normalizePath(file.path(images.info$path, 
                                                             basename(images.info$file))), 
                              to = file.path(images.info$rename.path, 
                                             images.info$rename.file),
                              ...)
    # Add status column to image into
    images.info$status <- status
    
    # Check
    if(any(!status)) {
      warning("rename_mda: At least some files were not renamed, see warnings!")
    } else {
      cat(paste("rename_mda: It seems the renaming went well :) check your output directory at:", 
                unique(images.info$rename.path), "\n"))
    }
  } else {
    cat(paste("rename_mda: rename.function set to NULL; no images were renamed to path:", 
              unique(images.info$rename.path), "\n"))
  }
  
  # Return mapping.
  return(list(
    rename.path=unique(images.info$rename.path),
    status=images.info
  ))
}
