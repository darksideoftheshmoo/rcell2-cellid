#' A function to copy the spos XLSX spreadsheet template to the current working directory.
#' 
#' The template file is included in \code{cell2.cellid}.
#' 
#' @param file_name File name for the downloaded workflow template.
#' 
#' @export
 get_spos_template <- function(file_name = "stage_positions.xlsx"){
  
  if(file.exists(file_name)) stop("file '", file_name, "' exists.")
  
  workflow.file <- system.file(
    "stage_positions.xlsx",
    package = "rcell2.cellid"
  )
  
  if(file.exists(workflow.file)){
    file.copy(from = workflow.file, to = file_name)
    cat(paste0("Spreadsheet template copied to the current working directory at: '", getwd(), file_name, "'"))
  } else {
    stop("file", file_name, "not found in package. Please report.")
  }
}


#' Generate stage position data for a 384 well-plate experiment
#' 
#' @details
#' This helper function generates a ".STG" file that can be imported into MetaMorph's MDA, to load stage positions for an experiment.
#' 
#' The function requires an XLSX spreadsheet as an input, which must contain sheets with:
#' 
#' 1. "pdata" with the position metadata,
#' 
#' 2. "well_order" with a representation of the plate indicating the order in which to take pictures,
#' 
#' 3. and "well_images" with a representation of the plate indicating how many images to acquire in each well.
#' 
#' A template XLSX file can be obtained with \code{get_spos_template}.
#' 
#' @param spreadsheet_path Path to the an XLSX spreadsheet (e.g. "stage_positions.xlsx") with the following sheets: "pdata", "well_order", "well_images." A template file can be obtained with \code{get_spos_template}.
#' @param spreadsheet_path Path where the STG file will be saved. This file is meant to be imported with MetaMorph's MDA.
#' @param well_sep Center-to-center distance between wells in the plate (in microns).
#' @param well_width Width of the bottom of the wells (in microns).
#' @param fov_width Width of the field of view of the final image (in microns).
#' @param origin_at_corner Corner of the well at which the coordinate origin was set.
#' @param origin_at_pos Stage position that contains the origin (at the corner selected by \code{origin_at_corner}).
#' @param af_offset_default Default value for the hardware auto-focus offset.
#' @param z_default Default value for the Z-coordinate (i.e. objective lens height).
#' 
#' @export
#' @import readxl ggplot2 dplyr tidyr
make_stage_list <- function(
  spreadsheet_path,
  stg_output_path="stage_list.STG",
  well_sep = 4500,
  well_width = 3300,
  fov_width=500,
  origin_at_corner = "top-right",
  origin_at_pos = 1,
  af_offset_default = 0,
  z_default = 0){
  
  pdata <- readxl::read_xlsx(spreadsheet_path, "pdata")
  
  well_order <- readxl::read_xlsx(spreadsheet_path, "well_order")
  names(well_order)[1] <- "row"
  well_order <- well_order |> 
    pivot_longer(-row, names_to = "col", values_to = "order") |> 
    filter(order > 0) |> 
    arrange(order)
  
  if(!all(diff(well_order$order) == 1)) stop("Ordering is not consecutive.")
  
  well_images <- readxl::read_xlsx(spreadsheet_path, "well_images")
  names(well_images)[1] <- "row"
  well_images <- well_images |> 
    pivot_longer(-row, names_to = "col", values_to = "images") |> 
    filter(images > 0)
  
  wells <- 
    left_join(well_order, well_images, 
              by = join_by(row, col)) |> 
    arrange(order) 
  
  # positions <- wells |> 
  #   mutate(pos_count = cumsum(images)) |> 
  #   mutate(start_pos = pos_count - images + 1,
  #          end_pos = start_pos + images - 1)
  
  positions <- wells |> 
    group_by(row, col, order) |> 
    reframe(fov=1:max(images)) |> 
    arrange(order, fov) |> 
    mutate(pos = 1:n())
  
  if(!all(positions$pos == pdata$pos)) stop("Position indexes in pdata do not match the well images and ordering.")
  
  if(origin_at_corner != "top-right") stop("Only top-right origins supported ATM.")
  if(origin_at_pos != 1) stop("The origin can only be in the well of the first position ATM.")
  
  half_n_fov <- floor((well_width %/% fov_width)/2)
  
  stage_coords <-  positions |> 
    select(pos,order,row,col,fov) |> 
    dplyr::rename(well=order) |> 
    mutate(row_i = as.integer(factor(row, levels=LETTERS))) |> 
    mutate(col_i = as.integer(col)) |> 
    mutate(
      row_i = row_i - min(row_i),
      col_i = col_i - min(col_i),
      fov_i = fov - min(fov)
    ) |> 
    mutate(
      x = (-well_width/2) + (col_i * well_sep) + (fov_i %% half_n_fov) * fov_width,
      y = (-well_width/2) - (row_i * well_sep) + (fov_i %/% half_n_fov) * fov_width,
      z = z_default,
      af_offset = af_offset_default
    ) |> 
    group_by(row, col) |> 
    mutate(
      x = x - max(fov_i %% half_n_fov)/2 * fov_width,
      y = y - max(fov_i %/% half_n_fov)/2 * fov_width
    ) |> 
    ungroup()
  
  plt <- 
    stage_coords |> 
    mutate(
      well = as.factor(well),
      fov = as.factor(fov)
    ) |> 
    ggplot() +
    geom_point(aes(x,y),color="red", shape=23, size=5, data=data.frame(x=0,y=0)) +
    geom_point(aes(x,y),color="red", shape=19, size=2, data=data.frame(x=0,y=0)) +
    geom_path(aes(x,y)) +
    geom_point(aes(x,y,color=well), size=5) +
    scale_color_discrete() +
    ggtitle("Generated stage coordinates for each imaging position",
            "The red diamond indicates the location of the origin (" |> 
              paste0(origin_at_corner, " corner, at position ", origin_at_pos, ").")) +
    theme_minimal()
  
  stg_file_header <- paste(
    '"Stage Memory List", Version 6.0',
    '0, 0, 0, 0, 0, 0, 0, "um", "um"',
    '0',
    nrow(positions),
    sep = "\n"
  )
  stg_file_data <- stage_coords |> 
    arrange(pos) |> 
    select(pos, x, y, z, af_offset) |> 
    mutate(pos = paste0("Position", pos)) |> 
    mutate(z2 = z) |> 
    mutate(
      v1 = FALSE,
      v2 = -9999,
      v3 = TRUE,
      v4 = TRUE,
      v5 = 0,
      v5 = -1,
      v6 = ""
    )
  write(x = stg_file_header, file = stg_output_path, append = F)
  write.table(stg_file_data, sep = ", ", file=stg_output_path, col.names = F, append = T, row.names = F)
  
  cat(paste0("Wrote the following content to: ", stg_output_path, "\n"))
  readLines(stg_output_path) |> paste(collapse = "\n") |> cat()

  return(list(
    stg_output_path=stg_output_path,
    positions_plot=plt, 
    stage_coords=stage_coords
  ))
}
