#' Visualise locations of gene conversion events
#'
#' @param annotation data.frame, output from \code{batchConvertAnalysis}.
#' @param repertoire A named vector of strings storing nucleotide sequences observed from a repertoire. The \code{names} attribute of the vector stores the sequence identifiers.
#' @param thickness A numeric vector. Controls the thickness of lines plotted. The first number corresponds to the line representing the region from \code{start} to \code{end} (ie the narrow definition). The second number corresponds to lines drawn fro the broad definition of events including the padded sequences 5'/3' to the start/end points.
#' @param show.sequence.names Do you want the sequence identifiers to be printed on the vertical axis of the output plot? Recommended to set to \code{FALSE} if many sequences are supplied in \code{annotation}.
#' @param highlight.genes optional, a character vector consisting of gene names which are to be highlighted with different colours in the output plot.
#'
#' @description This function visualises gene conversion events along the V gene listed in \code{annotation}. The start-end ranges are represented as thick lines and the 'broader' definition including the padded 5'/3' sequence stretched are represented as thin lines. Vertical lines indicate boundaries between framework regions and CDR.
#'
#' @return A ggplot2 object. See description for details. If \code{highlight.genes} is not \code{NULL} and the given genes can be found in \code{annotation$gene}, gene conversion events using these genes will be highlighted in different colours.
#'
#' @import ggplot2
#' @export plotEvents
plotEvents <- function(annotation, repertoire, thickness = c(0.8, 0.3),
                       show.sequence.names = FALSE, highlight.genes = NULL)
{
  if( !is.logical( show.sequence.names ) )
    stop("'show.sequence.names' must be either TRUE or FALSE.")
  if( !is.numeric( thickness ) )
    stop("'thickness' should be a numeric vector.")
  if( length( thickness ) != 2 )
    stop("'thickness' should be a vector of 2 numbers. See documentation for details.")
  if( show.sequence.names ){
    show.sequence.names = element_text()
  } else {
    show.sequence.names = element_blank()
  }
  # get 'soft' start and end positions
  tb <- adjustBroadPos( annotation, repertoire )
  if( !is.null( highlight.genes ) ){
    if( !is.character( highlight.genes ) ) stop("'highlight.genes' should be a vector of gene names to be highlighted in the plot.")
    tb <- tidyr::separate_rows(tb, as.symbol("gene"), sep=";")
    tb$gene_included <- sapply( tb$gene, function(x){
      if(is.na(x)) return(NA)
      # if allele numbers (i.e. *01, *02 etc) in highlight.genes than do exact match
      # otherwise just match by the gene name
      if( any(grepl("*", highlight.genes, fixed = TRUE)) ){
        if(x %in% highlight.genes) return(x)
      } else {
        o <- c(stringi::stri_match(x, regex = paste0(highlight.genes, "\\*")))
        if(all(is.na(o))) return(NA) else return( o[!is.na(o)])
      }
      return(NA)
    })
  }
  if( all( is.na(tb$gene_included) ) ){
    warning("None of the genes in highlight.genes are found in 'annotation'.")
    highlight.genes <- NULL
  }
  # adjust the start/end points to fit between [1, 330] (i.e. AA pos 1 - 110,
  # roughly the V gene boundary) - just in case.
  tb$start <- replace(tb$start, which(tb$start < 1), 1)
  tb$soft_start <- replace(tb$soft_start, which(tb$soft_start < 1), 1)
  tb$end <- replace(tb$end, which(tb$end > 330), 330)
  tb$soft_end <- replace(tb$soft_end, which(tb$soft_end > 330), 330)
  g <- ggplot() + geom_vline(xintercept = c(26, 38, 55, 65, 104, 110) * 3,
                             linetype = "dashed", colour = "grey60")
  if( !is.null( highlight.genes ) ){
    g <- g + geom_segment(data = tb[is.na(tb$gene_included), ],
                          aes_string(y = "SeqID", yend = "SeqID", x = "soft_start", xend = "soft_end", colour = "gene_included"),
                          size = thickness[2]) +
      geom_segment(data = tb[is.na(tb$gene_included), ],
                   aes_string(y = "SeqID", yend = "SeqID", x = "start", xend = "end", colour = "gene_included"),
                   size = thickness[1]) +
      geom_segment(data = tb[!is.na(tb$gene_included), ],
                          aes_string(y = "SeqID", yend = "SeqID", x = "soft_start", xend = "soft_end", colour = "gene_included"),
                          size = thickness[2]) +
      geom_segment(data = tb[!is.na(tb$gene_included), ],
                   aes_string(y = "SeqID", yend = "SeqID", x = "start", xend = "end", colour = "gene_included"),
                   size = thickness[1]) +
      scale_colour_discrete(labels = c(levels(factor(highlight.genes)), "Others"),
                            name = "genes", na.value = "grey40")
  } else {
    g <- g + geom_segment(data = tb, aes_string(y = "SeqID", yend = "SeqID", x = "soft_start", xend = "soft_end"), size = thickness[2]) +
      geom_segment(data = tb, aes_string(y = "SeqID", yend = "SeqID", x = "start", xend = "end"), size = thickness[1])
  }
#  g <- g + geom_rect(data=data.frame(xmin=c(26, 55, 104) * 3 + 1,
#                                     xmax=c(38, 65, 110) * 3),
#                     aes_string(xmin="xmin",xmax="xmax",ymin=1,ymax=length(unique(tb$SeqID))),
#                     fill="grey70", alpha=0.5, inherit.aes = FALSE)
  g <- g + cowplot::theme_cowplot() + xlab("DNA position") + ylab("Sequences") +
    theme(axis.text.y = show.sequence.names, axis.ticks.y = element_blank()) +
    scale_x_continuous(breaks = c(1, 101, 201, 301), limits = c(0, 330))
  g
}

#' Histograms of gene conversion events
#'
#' @param tb data.frame, output from \code{batchConvertAnalysis}.
#' @param repertoire A named vector of strings storing nucleotide sequences observed from a repertoire. The \code{names} attribute of the vector stores the sequence identifiers.
#' @param column1 optional, character corresponding to a column name in \code{tb}. If supplied, calculation of coverage will be considered separately for every unique value in \code{column1} of the data frame \code{tb} so that histograms can be plotted separately.
#' @param column2 optional, character corresponding to a column name in \code{tb}. If supplied, calculation of coverage will be considered separately for every unique value in \code{column2} of the data frame \code{tb} so that histograms can be plotted separately.
#' @param type character, either 'coverage' (i.e. how often is a given DNA position included within gene conversion events?) or 'start' (i.e. how often do gene conversion events start at a given DNA position?).
#' @param definition character. Indicates which definition of gene conversion event boundaries are to be used. Can take the values of any one of 'narrow' (the strict start/end points given by \code{batchConvertAnalysis}), 'broad' (the start/end padded with nucleotides 5'/3' to the given start/end) or 'both' (both 'narrow' and 'end' will be considered and two sets of plots will be generated.)
#'
#' @description This function visualises histograms showing the coverage or starting point of gene conversion events in a given set of annotation generated from \code{batchConvertAnalysis}, using definition of events.
#'
#' @return A ggplot2 object showing the histograms, with x-axis showing the DNA position and y-axis showing the percentage of gene conversion events which cover/start at the given position. Framework (FW) or CDR positions are coloured differently. If \code{column1} and/or \code{column2} are provided separate histograms will be plotted for different data partitions specified along the two given columns present in \code{tb}, arranged side-by-side.
#'
#' @import ggplot2
#' @importFrom stats as.formula
#' @export plotCoverage
plotCoverage <- function(tb, repertoire, column1 = NULL, column2 = NULL,
                         type = c("coverage", "start"),
                         definition = c("narrow", "broad", "both"))
{
  if( ! type %in% c("coverage", "start") )
    stop("'type' must be one of 'coverage' or 'start'. See documentation for details.")
  if( ! definition %in% c("narrow", "broad", "both") )
    stop("'definition' must be one of 'narrow', 'broad', or 'both'. See documentation for details.")
  tb <- adjustBroadPos( tb, repertoire )
  if( definition == 'both' ){
    coverage <- list(getCoverage(tb, column1, column2, type = type,
                                 start_col = "start", end_col = "end"),
                     getCoverage(tb, column1, column2, type = type,
                                 start_col = "soft_start", end_col = "soft_end"))
    coverage[[1]]$type <- "narrow event"
    coverage[[2]]$type <- "broad event"
    g <- list(
      ggplot(coverage[[1]], aes_string(x = "position", y = "value")) +
        geom_bar(aes_string(fill = "region"), stat = "identity") +
        scale_fill_manual(values = c("CDR" = "black", "FW" = "grey")) +
        cowplot::theme_cowplot() +
        scale_y_continuous(labels = scales::percent, name = "% gene conversion events") +
        scale_x_continuous(breaks = c(1, 101, 201, 301), name = "DNA position"),
      ggplot(coverage[[2]], aes_string(x = "position", y = "value")) +
        geom_bar(aes_string(fill = "region"), stat = "identity") +
        scale_fill_manual(values = c("CDR" = "black", "FW" = "grey")) +
        cowplot::theme_cowplot() +
        scale_y_continuous(labels = scales::percent, name = "% gene conversion events") +
        scale_x_continuous(breaks = c(1, 101, 201, 301), name = "DNA position")
    )
    g[[1]] <- g[[1]] + ggtitle("Narrow")
    g[[2]] <- g[[2]] + ggtitle("Broad")
    if( !any( is.null(c(column1, column2)) ) ){
      if( is.null(column1) ) column1 <- "."
      if( is.null(column2) ) column2 <- "."
      g[[1]] <- g[[1]] + facet_grid(as.formula(paste( column1, "~", column2 ) ))
      g[[2]] <- g[[2]] + facet_grid(as.formula(paste( column1, "~", column2 ) ))
    }
    return(
      cowplot::plot_grid(
        plotlist = g, ncol = 1, align = "v", axis = "lr"
      )
    )
  } else if (definition == 'narrow') {
    coverage <- getCoverage(tb, column1, column2, type = type,
                            start_col = "start", end_col = "end")
  } else if (definition == 'broad') {
    coverage <- getCoverage(tb, column1, column2, type = type,
                            start_col = "soft_start", end_col = "soft_end")
  }
  g <- ggplot(coverage, aes_string(x = "position", y = "value")) +
    geom_bar(aes_string(fill = "region"), stat = "identity") +
    scale_fill_manual(values = c("CDR" = "black", "FW" = "grey")) +
    cowplot::theme_cowplot() +
    scale_y_continuous(labels = scales::percent, name = "% gene conversion events") +
    scale_x_continuous(breaks = c(1, 101, 201, 301), name = "DNA position")
  if( !any( is.null(c(column1, column2)) ) ){
    if( is.null(column1) ) column1 <- "."
    if( is.null(column2) ) column2 <- "."
    g <- g + facet_grid(as.formula(paste( column1, "~", column2 ) ))
  }
  return( g )
}
