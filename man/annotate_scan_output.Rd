% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cellid_parameter_scan.R
\name{annotate_scan_output}
\alias{annotate_scan_output}
\title{Annotate parameter scan output for reviewing the results in ImageJ}
\usage{
annotate_scan_output(
  scan.results,
  annotate.channels = "BF.out",
  preserve_source_imgs = TRUE,
  in.place = FALSE,
  annotated.imgs.dir = "annotated",
  annotation.font = NULL
)
}
\arguments{
\item{scan.results}{The full result from \code{parameter_scan}.}

\item{annotate.channels}{Vector of channels which should be annotated.}

\item{preserve_source_imgs}{If TRUE, an error will be raised when non-output (segmented) images are found in the input.}

\item{in.place}{If TRUE, the provided images will be annotated and overwritten. Else, the annotated images will be saved to a subdirectory of \code{test.dir} named by \code{annotated.imgs.dir}.}

\item{annotated.imgs.dir}{If \code{in.place} is FALSE, the annotated images will be saved to a sub-directory of \code{test.dir} named by \code{annotated.imgs.dir}, and suffixed by the channel names.}

\item{annotation.font}{Font for the annotations. A mono-spaced font is recommended.}
}
\value{
ImageJ macros to open the annotated images as a virtual stack.
}
\description{
Warning: images in the scan result will be annotated and overwritten. Be sure to double check that your original "source" images are _not_ in the data frame (i.e. pass only output ".out.tif" images to this function).
}
