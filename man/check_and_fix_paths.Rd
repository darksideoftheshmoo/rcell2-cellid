% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cellid_wrapper.R
\name{check_and_fix_paths}
\alias{check_and_fix_paths}
\title{Check and fix image paths}
\usage{
check_and_fix_paths(images, path = NULL)
}
\arguments{
\item{images}{Cell-ID data-frame with image paths, as loaded by \code{get_cell_data}.}

\item{path}{Vector of paths to alternate data directories, where the images are actually stored. If left NULL, the function will raise an error when images are missing.}
}
\description{
This function tries to find files for images at the paths found in the provided data frame, and optionally fixes them if an alternative path is provided.
}
