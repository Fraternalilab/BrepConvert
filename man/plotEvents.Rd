% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{plotEvents}
\alias{plotEvents}
\title{Visualise locations of gene conversion events}
\usage{
plotEvents(
  annotation,
  repertoire,
  thickness = c(0.8, 0.3),
  show.sequence.names = FALSE,
  highlight.genes = NULL
)
}
\arguments{
\item{annotation}{data.frame, output from \code{batchConvertAnalysis}.}

\item{repertoire}{A named vector of strings storing nucleotide sequences observed from a repertoire. The \code{names} attribute of the vector stores the sequence identifiers.}

\item{thickness}{A numeric vector. Controls the thickness of lines plotted. The first number corresponds to the line representing the region from \code{start} to \code{end} (ie the narrow definition). The second number corresponds to lines drawn fro the broad definition of events including the padded sequences 5'/3' to the start/end points.}

\item{show.sequence.names}{Do you want the sequence identifiers to be printed on the vertical axis of the output plot? Recommended to set to \code{FALSE} if many sequences are supplied in \code{annotation}.}

\item{highlight.genes}{optional, a character vector consisting of gene names which are to be highlighted with different colours in the output plot.}
}
\value{
A ggplot2 object. See description for details. If \code{highlight.genes} is not \code{NULL} and the given genes can be found in \code{annotation$gene}, gene conversion events using these genes will be highlighted in different colours.
}
\description{
This function visualises gene conversion events along the V gene listed in \code{annotation}. The start-end ranges are represented as thick lines and the 'broader' definition including the padded 5'/3' sequence stretched are represented as thin lines. Vertical lines indicate boundaries between framework regions and CDR.
}
