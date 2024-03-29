% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cellid_wrapper.R
\name{arguments_check}
\alias{arguments_check}
\title{Check if output files have already been created}
\usage{
arguments_check(arguments_df, check_fail = F)
}
\arguments{
\item{arguments_df}{The table generated by \code{arguments}.}

\item{check_fail}{Whether to produce an error if any output file or directory already exists.}
}
\value{
TRUE if no existing output files were found.
}
\description{
Check if output files have already been created
}
\details{
Emits warnings, and optionally an error, to prevent accidental overwriting of existing output files.

Also checks if Cell-IDs t.frame numbering will match the t.frame in the file's names.

Also checks that the channel identifier is 3 characters long, as required by Cell-ID to distinguish imaging channels.
}
