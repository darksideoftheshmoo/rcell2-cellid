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
  
  if(!requireNamespace("XML")){
    warning("tiff_plane_info requires xmlToList from the 'XML' package, which is not installed. Aborting.")
    return(NULL)
  }
  
  if(!requireNamespace("ijtiff")){
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
        "tiff_plane_info: error while parsing metadata from file '",
        basename(path), "'. Returning NULL."))
      return(NULL)
    })
  
  return(result)
}

#' Plot positions and field of view overlaps
#' 
#' @param image_list The output from the \code{arguments} function. Consider fltering to include only the first position.
#' @param magnification Total image magnification.
#' @param ccd_pixel_size_microns Physical size of the camera's pixel.
#' @param well_size_microns_x Physical size of the well in the plate (used to draw breaks in the plot).
#' @param well_size_microns_y Physical size of the well in the plate (used to draw breaks in the plot).
#' @param image_width_x Size of the images in pixels.
#' @param image_height_y Size of the images in pixels.
#' @import dplyr ggplot2
#' @export
#' @examples 
#' plot <- rcell2.cellid::arguments(metamorph_pics_dir, file.pattern = "^(BF|[A-Z]FP)_Position(\\d+)_time(\\d+).tif$") |> 
#'   filter(t.frame == min(t.frame)) |> 
#'   plot_pos_overlaps()
#' 
plot_pos_overlaps <- function(
    image_list,
    # file.pattern = "^(BF|[A-Z]FP)_Position(\\d+)_time(\\d+).tif$",
    magnification = 40,            # 40x
    ccd_pixel_size_microns = 6.45, # 6.45 um
    well_size_microns_x = 4500,    # 4500 um
    well_size_microns_y = 4500,    # 4500 um
    image_width_x = 1376,
    image_height_y = 1040
    ){
  
  # Path to images.
  metamorph_pics <- paste0(image_list$path, "/", image_list$bf) |> unique()
  
  # Get the metadata of all images:
  plane_info_df <- setNames(metamorph_pics, basename(metamorph_pics)) %>% 
    # Get metadata
    lapply(rcell2.cellid::tiff_plane_info) %>% 
    bind_rows(.id="image") %>% 
    # Fix variable names and types
    mutate(variable = make.names(id)) %>% 
    filter(grepl("stage.position", variable)) %>% 
    mutate(value = as.numeric(value)) %>% 
    # Make it wider
    pivot_wider(id_cols = "image", names_from = "variable", values_from = "value")
  
  # Generate "images" dataframe.
  images <- arguments_to_images(arguments = image_list)
  
  # Microscope details.
  # magnification <- 40 # 40x
  # ccd_pixel_size_microns <- 6.45 # 6.45 um
  
  # Well size
  # well_size_microns_x <- 9000 / 2
  # well_size_microns_y <- 9000 / 2
  
  # Calculate "Field Of View" size.
  fov_size_microns_x <- image_width_x * ccd_pixel_size_microns / magnification  # um in the "X/width" direction
  fov_size_microns_y <- image_height_y * ccd_pixel_size_microns / magnification  # um in the "Y/height" direction
  
  x_min <- min(plane_info_df$stage.position.x) - well_size_microns_x / 2
  x_max <- max(plane_info_df$stage.position.x) + well_size_microns_x / 2
  y_min <- min(plane_info_df$stage.position.y) - well_size_microns_x / 2
  y_max <- max(plane_info_df$stage.position.y) + well_size_microns_x / 2
  
  images_bf <- select(images, image, pos, t.frame, channel) |> filter(channel == "BF") |> unique()
  
  # Plot fields of view to check for overlaps visually.
  plt <- plane_info_df %>% 
    left_join(images_bf, by = "image") |> 
    arrange(pos) %>% 
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
  
  print(plt)
  
  return(invisible(list(
    plt=plt,
    plane_info_df=plane_info_df,
    images=images
  )))
}
