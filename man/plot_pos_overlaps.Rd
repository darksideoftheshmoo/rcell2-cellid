% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_tiff_info.R
\name{plot_pos_overlaps}
\alias{plot_pos_overlaps}
\title{Plot positions and field of view overlaps}
\usage{
plot_pos_overlaps(
  image_list = NULL,
  images = NULL,
  channels = "BF",
  magnification = 40,
  ccd_pixel_size_microns = 6.45,
  well_size_microns_x = 4500,
  well_size_microns_y = 4500,
  image_width_x = 1376,
  image_height_y = 1040,
  print_plot = T
)
}
\arguments{
\item{image_list}{The output from the \code{arguments} function. Consider filtering to include only the first position.}

\item{images}{The images dataframe as loaded from Cell-ID's output by \code{get_cell_data}.}

\item{channels}{Vector of channel IDs (e.g. "BF").}

\item{magnification}{Total image magnification.}

\item{ccd_pixel_size_microns}{Physical size of the camera's pixel.}

\item{well_size_microns_x}{Physical size of the well in the plate (used to draw breaks in the plot).}

\item{well_size_microns_y}{Physical size of the well in the plate (used to draw breaks in the plot).}

\item{image_width_x}{Size of the images in pixels.}

\item{image_height_y}{Size of the images in pixels.}

\item{print_plot}{Print the plot.}
}
\description{
Example:
plot <- rcell2.cellid::arguments(metamorph_pics_dir, file.pattern = "^(BF|[A-Z]FP)_Position(\\d+)_time(\\d+).tif$") |> 
  filter(t.frame == min(t.frame)) |> 
  plot_pos_overlaps()
}
