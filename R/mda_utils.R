#' A function to copy the spos XLSX spreadsheet template to the current working directory.
#' 
#' The template file is included in \code{cell2.cellid}.
#' 
#' @param save_to File name for the downloaded workflow template.
#' @param overwrite Flag to overwrite the existing file if found.
#' 
#' @export
 get_spos_template <- function(save_to = "stage_positions.xlsx", overwrite = FALSE){
  
  if(file.exists(save_to) & !overwrite) stop("file '", save_to, "' already exists.")
  
  workflow.file <- system.file(
    "stage_positions.xlsx",
    package = "rcell2.cellid"
  )
  
  if(file.exists(workflow.file)){
    cat(paste0("Spreadsheet template will be saved to: '", file.path(getwd(), save_to), "'"))
    # Raise an error if it fails.
    stopifnot(
      file.copy(from = workflow.file, to = save_to, overwrite = overwrite)
    )
  } else {
    stop("File", save_to, "not found in package. Please report.")
  }
  
  return(save_to)
}


#' Parse a grid sheet into long format
grid_to_long <- function(sheet_data, rows_to="row", names_to="col", values_to="row"){
  names(sheet_data)[1] <- rows_to
  sheet_data <- sheet_data |> 
    pivot_longer(-all_of(rows_to), names_to=names_to, values_to=values_to)
  return(sheet_data)
}

#' Extract well order and image numbers from the template spreadsheet
make_wells <- function(spreadsheet_path,
                       well_order_sheet="well_order", 
                       well_images_sheet="well_images"){
  
  well_order <- readxl::read_xlsx(spreadsheet_path, well_order_sheet) |> 
    grid_to_long(names_to = "col", values_to = "order") |> 
    filter(order > 0) |> 
    arrange(order)
  
  if(!all(diff(well_order$order) == 1)) stop("Ordering is not consecutive.")
  
  well_images <- readxl::read_xlsx(spreadsheet_path, well_images_sheet) |> 
    grid_to_long(names_to = "col", values_to = "images") |> 
    filter(images > 0)
  
  wells <- 
    left_join(well_order, well_images, 
              by = join_by(row, col)) |> 
    arrange(order) 
  
  return(wells)
}

#' Prepare a pdata dataframe from a template spreadsheet
#' 
#' Extract metadata information from a template spreadsheet (e.g. as the one provided by \code{get_spos_template}), joining it to the "pdata" table in that spreadsheet.
#' 
#' @inheritParams make_stage_list
#' @param pdata_sheets Extract additional metadata from these sheets in the template spreadsheet file.
#' @param plot_metadata If TRUE, plot the extracted metadata.
#' @param write_pdata If TRUE, write a "pdata.csv" file next to the original template file.
#' @export
#' @import ggplot2 dplyr tidyr readxl
make_pdata <- function(spreadsheet_path, pdata_sheets=c(), plot_metadata=TRUE, write_pdata=TRUE){
  
  pdata <- readxl::read_xlsx(spreadsheet_path, "pdata")
  
  positions <- spreadsheet_path |> 
    make_wells() |> 
    make_positions()
  
  if(!all(sort(positions$pos) == sort(pdata$pos))) stop("Position indexes in pdata do not match the well images and ordering.")
  
  pdata <- pdata |> 
    left_join(positions, by = "pos")
  
  for(sheet in pdata_sheets){
    metadata <- readxl::read_xlsx(spreadsheet_path, sheet) |> 
      grid_to_long(values_to=sheet)
    
    pdata <- pdata |> 
      left_join(metadata, by=c("row", "col"))
  }
  
  pdata <- pdata |> 
    mutate_at(.vars = c("row", "col"), .funs = as.ordered)
  
  if(plot_metadata){
    for(sheet in pdata_sheets){
      plt <- pdata |> 
        ggplot() +
        geom_tile(aes(col, row, fill=as.ordered(.data[[sheet]]))) + 
        scale_y_discrete(limits=rev) +
        labs(fill=sheet) + theme_minimal()
      
      print(plt)
    }
  }
  
  if(isTRUE(write_pdata)){
    write.csv(x = pdata, 
              file = file.path(dirname(spreadsheet_path), "pdata.csv"))
  }
  
  return(pdata)
}

#' Add 'pos' and 'fov' fields to the wells data frame
make_positions <- function(wells){
  
  positions <- wells |> 
    group_by(row, col, order) |> 
    reframe(fov=1:max(images)) |> 
    arrange(order, fov) |> 
    mutate(pos = 1:n())
  
  return(positions)
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
#' @param stg_output_path Path where the STG file will be saved. This file is meant to be imported with MetaMorph's MDA.
#' @param calib_stg_path Path where the calibration STG can be found saved. NULL by default (disabled).
#' @param well_sep Center-to-center distance between wells in the plate (in microns).
#' @param well_width Width of the bottom of the wells (in microns).
#' @param fov_width Width of the field of view of the final image (in microns).
#' @param origin_at_corner Corner of the well at which the coordinate origin was set.
#' @param af_offset_default Default value for the hardware auto-focus offset.
#' @param z_default Default value for the Z-coordinate (i.e. objective lens height).
#' 
#' @export
#' @import readxl ggplot2 dplyr tidyr
make_stage_list <- function(
  spreadsheet_path,
  stg_output_path="stage_list.STG",
  print_output=FALSE,
  calib_stg_path=NULL,
  well_sep = 4500,
  well_width = 3300,
  fov_width=500,
  origin_at_corner = "top-right",
  # origin_at_pos = 1,
  # @param origin_at_pos Stage position that contains the origin (at the corner selected by \code{origin_at_corner}).
  af_offset_default = 0,
  z_default = 0){
  
  pdata <- readxl::read_xlsx(spreadsheet_path, "pdata")
  
  wells <- make_wells(spreadsheet_path)
  
  positions <- make_positions(wells)
  
  if(!all(sort(positions$pos) == sort(pdata$pos))) stop("Position indexes in pdata do not match the well images and ordering.")
  
  if(origin_at_corner != "top-right") stop("Only top-right origins supported ATM.")
  # if(origin_at_pos != 1) stop("The origin can only be in the well of the first position ATM.")
  
  half_n_fov <- floor((well_width %/% fov_width)/2)
  
  stage_idxs <-  positions |> 
    select(pos,order,row,col,fov) |> 
    dplyr::rename(well=order) |> 
    mutate(row_i = as.integer(factor(row, levels=LETTERS))) |> 
    mutate(col_i = as.integer(col)) |> 
    mutate(
      row_i = row_i - min(row_i),
      col_i = col_i - min(col_i),
      fov_i = fov - min(fov)
    )
  
  origin_at <- stage_idxs |> 
    filter(row_i == min(row_i)) |> 
    filter(col_i == min(col_i)) |> 
    with(list(pos=pos, row_i=row_i, col_i=col_i))
  
  origin_at_pos <- origin_at$pos
  
  stage_coords <- stage_idxs |> 
    mutate(
      row_i = row_i - origin_at$row_i,
      col_i = col_i - origin_at$col_i
    ) |> 
    mutate(
      # The stage coordinates use more negative numbers to move towards the right of a well-plate.
      x = (well_width/2) - (col_i * well_sep) - (fov_i %% half_n_fov) * fov_width,
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
  
  # Adjust with calibration
  if(!is.null(calib_stg_path)){
    # Get the default calibration file.
    if(isTRUE(calib_stg_path)){
      calib_stg_path <- system.file("stage_positions-calib_tilt_B2_top_left.STG", package = "rcell2.cellid")
    }
    # Adjust the Z coordinate.
    warning(paste("Adjusting stage coords with calibration file:", calib_stg_path))
    stage_coords <- adjust_stg_z(stage_coords, calib_stg_path = calib_stg_path, plot_error = F)
  }
  
  # Prepare plot objects.
  plt <- plot_stg_coords(stage_coords, origin_at_corner, origin_at_pos)
  plt_z <- plot_stg_coords_z(stage_coords, origin_at_corner, origin_at_pos)
  
  # Write position data in STG format
  stg_file_data <- write_stg(stage_coords, positions, stg_output_path)
  # Log.
  cat(paste0("Wrote the STG file to: ", normalizePath(stg_output_path), "\n"))
  if(print_output){
    readLines(stg_output_path) |> paste(collapse = "\n") |> cat()
  }

  # Done.
  return(list(
    stg_output_path=stg_output_path,
    positions_plot=plt, 
    xyz_plot=plt_z,
    stage_coords=stage_coords,
    stg_file_data=stg_file_data
  ))
}

plot_stg_coords <- function(stage_coords, origin_at_corner, origin_at_pos){
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
    geom_text(aes(x,y,label=pos), size=3) +
    # For the plot to render as a well-plate seen from above (as usual),
    # the X coordinates of the stage must be flipped.
    scale_x_reverse() +
    scale_color_discrete() +
    guides(colour="none") +
    ggtitle("Generated stage coordinates for each imaging position",
            "The red diamond indicates the location of the origin (" |> 
              paste0(origin_at_corner, " corner of position ", origin_at_pos, ").")) +
    theme_minimal()
  
  return(plt)
}

plot_stg_coords_z <- function(stage_coords, origin_at_corner, origin_at_pos){
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
    geom_point(aes(x,y,color=z), size=5) +
    # For the plot to render as a well-plate seen from above (as usual),
    # the X coordinates of the stage must be flipped.
    scale_x_reverse() +
    scale_color_viridis_c() +
    # guides(colour="none") +
    ggtitle("Generated stage coordinates for each imaging position",
            "The red diamond indicates the location of the origin (" |> 
              paste0(origin_at_corner, " corner of position ", origin_at_pos, ").")) +
    theme_minimal()
  
  return(plt)
}

#' Prepare stage data in STG format
make_stage_coords <- function(stage_coords){
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
  
  return(stg_file_data)
}

#' Write position data in STG format
#' 
#' @param stage_coords Stage coordinates data frame, as produced by \code{make_stage_list}.
#' @param positions Positions table as generated by \code{make_positions}.
#' @param stg_output_path Path where the STG file will be saved. This file is meant to be imported with MetaMorph's MDA.
#'
write_stg <- function(stage_coords, positions, stg_output_path){
  
  stg_file_header <- paste(
    '"Stage Memory List", Version 6.0',
    '0, 0, 0, 0, 0, 0, 0, "um", "um"',
    '0',
    nrow(positions),
    sep = "\n"
  )
  
  # Make output data frame.
  stg_file_data <- make_stage_coords(stage_coords)
  
  # Write the output.
  write(x = stg_file_header, file = stg_output_path, append = F)
  write.table(stg_file_data, sep = ", ", file=stg_output_path, col.names = F, append = T, row.names = F)
  
  return(stg_file_data)
}

#' Read STG file
#' 
#' @import ggplot2
#' @export
read_stg <- function(stg_list_path, plot_stage=TRUE){
  stg_list <- stg_list_path |> 
    read.csv(skip = 4, header = F)
  
  stg_list <- stg_list[1:5]
  
  names(stg_list) <- c(
    "position_name",
    "x", "y", "z",
    "af_offset"
  )
  
  if(plot_stage){
    plt <- stg_list |> 
      ggplot() +
      geom_path(aes(x,y)) +
      geom_point(aes(x,y,color=z), size=10, alpha =.5) +
      scale_x_reverse() + 
      scale_color_viridis_c() + theme_minimal()
    print(plt)
  }
  
  return(stg_list)
}

#' Fit a linear model to the XY coordinates and predict Z
model_stg_z <- function(calib_stg_list){
  m <- lm(z~x+y, data=calib_stg_list)
  return(m)
}

#' Use a calibration STG file to adjust the Z coordinates of a "stage_coords" data frame
#' @param stage_coords Stage coordinates data frame, as produced by "make_stage_list", or STG data as loaded by \code{read_stg}.
#' @param calib_stg_path Path to an STG file with calibration points. Their origin must be the same as the origin used in \code{stage_coords}.
#' @inheritParams make_stage_list
#' @export
adjust_stg_z <- function(
    stage_coords,
    calib_stg_path = system.file("stage_positions-calib_tilt_B2_top_left.STG", package = "rcell2.cellid"),
    origin_at_corner = "top-right",
    plot_error = FALSE
){
  if(origin_at_corner != "top-right") stop("Only top-right origins supported ATM.")
  
  calib_stg_list <- read_stg(calib_stg_path, plot_stage = F)
  m <- model_stg_z(calib_stg_list)
  xy <- stage_coords[,c("x", "y")]
  new_z <- predict(m, newdata=xy)
  stage_coords$z_err <- new_z - stage_coords$z
  stage_coords$z <- new_z
  
  if(plot_error){
    stage_coords |> adjust_stg_z_plot() |> print()
  }
  
  return(stage_coords)
}

adjust_stg_z_plot <- function(stage_coords){
  plt <- stage_coords |> 
    ggplot() +
    geom_path(aes(x,y)) +
    geom_point(aes(x,y,color=z_err), size=10, alpha =.5) +
    scale_x_reverse() + 
    scale_color_viridis_c() + theme_minimal()
  return(plt)
}

