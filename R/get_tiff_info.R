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
