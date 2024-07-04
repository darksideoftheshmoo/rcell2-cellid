#' Get "PlaneInfo" from a TIFF's tags (metadata)
#' 
#' @description 
#' 
#' This function parses the XML metadata in Metamprph's TIFF image files,
#' and retrieves the "plane info" elements, containing information on
#' position, exposure, etc.
#' 
#' Only the first frame's information will be shown. Use this function once per frame.
#' 
#' This function requires:
#' 
#' * xmlToList from the XML package.
#' 
#' * read_tags from the ijtiff package.
#' 
#' @inheritParams ijtiff::read_tags
#' @inheritDotParams XML::xmlToList
#' 
#' @import dplyr
#'
#' @return A dataframe with "prop" values of the TIFF's PlaneInfo metadata.
#' @export
#'
tiff_plane_info <- function(path, frames = 1, ...) {
  # @importFrom XML xmlToList
  # @importFrom ijtiff read_tags
  # path <- "~/Data/2023-06-08-screen_nuevas_cepas_Far1_3xmNG/data/images/renamed/YFP_Position1_time01.tif"
  
  if (length(frames) > 1) {
    warning("tiff_plane_info: only the first frame's information will be shown. Use this function once per frame.")
    frames <- frames[1]
  }
  
  if(!requireNamespace("XML", quietly = T)){
    warning("tiff_plane_info requires xmlToList from the 'XML' package, which is not installed. Aborting.")
    return(NULL)
  }
  
  if(!requireNamespace("ijtiff", quietly = T)){
    warning("tiff_plane_info requires read_tags from the 'ijtiff' package, which is not installed. Aborting.")
    return(NULL)
  }
  
  pic.tags <- ijtiff::read_tags(path = path, frames = frames)
  
  # Get metamorph's XML
  pic.description <- pic.tags$frame1$description
  
  if(is.null(pic.description)){
    warning("tiff_plane_info: Metamorph's plane-info XML metadata was not found in the TIFF's tags.")
    return(NULL)
  }
  
  # Catch XML errors, return NULL if it fails
  result <- tryCatch(
    expr = {
      # Parse XML to an R list
      description.xml <- XML::xmlToList(pic.description)
      # Convert plane info to a dataframe
      plane_info <- description.xml$PlaneInfo |>
        lapply(function(prop){
          # Allow "old and new" metadata.
          if(is.list(prop))
            return(prop$.attrs)  # Old (list and ".attrs")
          else
            return(prop) # New (vector)
        }) |>
        dplyr::bind_rows(.id = "prop_type")
      # Return
      plane_info
    },
    error = function(e){
      warning(e)
      warning(paste0(
        "\ntiff_plane_info: error while parsing metadata from file '",
        basename(path), "'. Returning NULL.\n"))
      return(NULL)
    })
  
  return(result)
}

#' Plot positions and field of view overlaps
#' 
#' Example:
#' plot <- rcell2.cellid::arguments(metamorph_pics_dir, file.pattern = "^(BF|[A-Z]FP)_Position(\\d+)_time(\\d+).tif$") |> 
#'   filter(t.frame == min(t.frame)) |> 
#'   plot_pos_overlaps()
#'   
#' @param image_list The output from the \code{arguments} function. Consider filtering to include only the first position.
#' @param images The images dataframe as loaded from Cell-ID's output by \code{get_cell_data}.
#' @param channels Vector of channel IDs (e.g. "BF").
#' @param magnification Total image magnification.
#' @param ccd_pixel_size_microns Physical size of the camera's pixel.
#' @param well_size_microns_x Physical size of the well in the plate (used to draw breaks in the plot).
#' @param well_size_microns_y Physical size of the well in the plate (used to draw breaks in the plot).
#' @param image_width_x Size of the images in pixels.
#' @param image_height_y Size of the images in pixels.
#' @param print_plot Print the plot.
#' @import dplyr ggplot2 tidyr
#' @export
plot_pos_overlaps <- function(
    image_list=NULL,
    images=NULL,
    channels = "BF",
    # file.pattern = "^(BF|[A-Z]FP)_Position(\\d+)_time(\\d+).tif$",
    magnification = 40,            # 40x
    ccd_pixel_size_microns = 6.45, # 6.45 um
    well_size_microns_x = 4500,    # 4500 um
    well_size_microns_y = 4500,    # 4500 um
    image_width_x = 1376,
    image_height_y = 1040,
    print_plot=T
    ){
  
  # Prepare the "images" dataframe.
  if (!is.null(image_list)) {
    # Convert Cell-ID arguments to images.
    available_channels <- image_list |> 
      with(ch) |> 
      unique()
    images <- rcell2.cellid::arguments_to_images(arguments = image_list) |> 
      dplyr::filter(channel %in% channels)
  } else if(is.null(images)){
    stop("Either 'image_list' or 'images' must be supplied.")
  } else {
    # Use images from Cell-ID's output.
    available_channels <- images |> 
      with(channels) |> 
      unique()
    images <- images |> 
      dplyr::filter(channel %in% channels)
  }
  
  if(!channels %in% available_channels) stop(paste(
    "Channels", setdiff(channels, available_channels), "are not available.",
    "Use any of the available channels instead:", available_channels
  ))
  
  # Path to images.
  metamorph_pics <- images |> with(file) |> unique()
  
  # Get the metadata of all images:
  plane_info_df <- setNames(metamorph_pics, basename(metamorph_pics)) %>% 
    # Get metadata
    lapply(rcell2.cellid::tiff_plane_info) %>% 
    dplyr::bind_rows(.id="image") %>% 
    # Fix variable names and types
    dplyr::mutate(variable = make.names(id)) %>% 
    # Get the stage positions.
    dplyr::filter(grepl("stage.position", variable)) %>% 
    dplyr::mutate(value = as.numeric(value)) %>% 
    # Make it wider
    tidyr::pivot_wider(id_cols = "image", names_from = "variable", values_from = "value")
  
  # Calculate "Field Of View" size.
  fov_size_microns_x <- image_width_x * ccd_pixel_size_microns / magnification  # um in the "X/width" direction
  fov_size_microns_y <- image_height_y * ccd_pixel_size_microns / magnification  # um in the "Y/height" direction
  
  x_min <- min(plane_info_df$stage.position.x) - well_size_microns_x / 2
  x_max <- max(plane_info_df$stage.position.x) + well_size_microns_x / 2
  y_min <- min(plane_info_df$stage.position.y) - well_size_microns_x / 2
  y_max <- max(plane_info_df$stage.position.y) + well_size_microns_x / 2
  
  images_bf <- select(images, image, pos, t.frame, channel) |> filter(channel %in% channels) |> unique()
  
  # Plot fields of view to check for overlaps visually.
  plt <- plane_info_df %>% 
    dplyr::left_join(images_bf, by = "image") |> 
    dplyr::arrange(pos) %>% 
    ggplot(aes(stage.position.x, stage.position.y, label = pos)) +
      geom_path(aes(group = t.frame)) +
      geom_rect(aes(xmin=stage.position.x-fov_size_microns_x/2, 
                    xmax=stage.position.x+fov_size_microns_x/2, 
                    ymin=stage.position.y-fov_size_microns_y/2, 
                    ymax=stage.position.y+fov_size_microns_y/2,
                    group = pos, fill = factor(pos)), alpha =.5)+
      geom_text(size=10) +
      facet_wrap(t.frame~channel) + guides(fill = "none") +
      
      scale_x_continuous(trans = "reverse", limits = c(x_max,x_min),
                         breaks=seq(0,-well_size_microns_x*10,by=-well_size_microns_x),
                         minor_breaks=NULL) + 
      scale_y_continuous(limits=c(y_min,y_max),
                         breaks=seq(0,-well_size_microns_y*10,by=-well_size_microns_y),
                         minor_breaks=NULL) +
      
      ggtitle("Physical stage coordinates v.s. Position index",
              "Compare the index numbers with the expected physical distrubution in the well plate.\nThe shaded areas around index numbers shuould not overlap with each other.")
  
  if(print_plot) print(plt)
  
  return(invisible(list(
    plt=plt,
    plane_info_df=plane_info_df,
    images=images
  )))
}


#' Load acquisition time information from MetaMorph's TIFF files
#' 
#' Time information is loaded from the XML inserted by MetaMorph into the TIFF files' metadata.
#' 
#' Time is loaded as POSIXlt and POSIXct into different columns (the 'ct' suffix indicates POSIXct). See \code{DateTimeClasses} for details.
#' 
#' @inheritParams plot_pos_overlaps
#' @param join_to_images If TRUE, the position/frame metadata form the images datafram will be joined to the extracted acquisition times.
#' @import dplyr tidyr
#' @export
get_time_info <- function(
    image_list=NULL,
    images=NULL,
    channels = "BF",
    join_to_images=T){
  
  # Prepare the "images" dataframe.
  if (!is.null(image_list)) {
    # Convert Cell-ID arguments to images.
    images <- rcell2.cellid::arguments_to_images(arguments = image_list) |> 
      dplyr::filter(channel %in% channels)
  } else if(is.null(images)){
    stop("Either 'image_list' or 'images' must be supplied.")
  } else {
    # Use images from Cell-ID's output.
    images <- images |> 
      dplyr::filter(channel %in% channels)
  }
  
  # Path to images.
  metamorph_pics <- images |> 
    with(file) |> unique()
  
  # Get the metadata of all images:
  tiff_info_df <- setNames(metamorph_pics, basename(metamorph_pics)) %>% 
    # Get metadata
    lapply(rcell2.cellid::tiff_plane_info) %>% 
    dplyr::bind_rows(.id="image") %>% 
    # Fix variable names and types
    dplyr::mutate(variable = make.names(id))
  
  time_info_df <- tiff_info_df %>% 
    # Get the stage positions.
    dplyr::filter(grepl("acquisition.time.local|modification.time.local", variable)) %>% 
    # 20231010 16:26:37.071
    dplyr::mutate(value = strptime(value, "%Y%m%d %H:%M:%OS")) %>% 
    # Make it wider
    tidyr::pivot_wider(id_cols = "image", names_from = "variable", values_from = "value")
  
  # Create columns with "POSIXct" time (required by ggplot).
  time_info_df <- time_info_df |> 
    mutate(acquisition.time.ct = as.POSIXct(acquisition.time.local),
           modification.time.ct = as.POSIXct(modification.time.local))
  
  # Join position and frame information.
  if(join_to_images){
    time_info_df <- time_info_df |> 
      dplyr::left_join(
        images |> dplyr::select(pos, t.frame, image) |> unique(),
        by = dplyr::join_by(image)
      )
  }
  
  return(time_info_df)
}