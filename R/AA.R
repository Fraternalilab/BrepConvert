#' Give 'broad' start/end points of gene conversion events
#'
#' @param tb data.frame, output from \code{batchConvertAnalysis}.
#' @param repertoire A named vector of strings storing nucleotide sequences observed from a repertoire. The \code{names} attribute of the vector stores the sequence identifiers.
#'
#' @description This function gives the 'broad' definition of gene conversion event boundaries, i.e. the start/end coordinates given in \code{batchConvertAnalysis} with added nucleotides corresponding to sequence stretches 5'/3' of start/end points identical between the observed repertoire sequence and the pseudogene. The output respects the same numbering scheme, respecting the imposed gaps ('.') in the input \code{repertoire}.
#'
#' @return A data.frame with these two columns added to \code{tb}:
#' \describe{
#'   \item{soft_start}{the position given in the \code{start} column, padded with the number of nucleotides given in \code{fiveprime_identical_length}.}
#'   \item{soft_end}{the position given in the \code{end} column, padded with the number of nucleotides given in \code{threeprime_identical_length}.}
#' }
#' @export adjustBroadPos
adjustBroadPos <- function(tb, repertoire)
{
  # adjust the soft start/end positions to give soft-starts/ends which
  # respect the gapped numbering
  tb$threeprime_identical_length <- as.numeric(tb$threeprime_identical_length)
  tb$fiveprime_identical_length <- as.numeric(tb$fiveprime_identical_length)
  # get ungapped coordinates of events for parsing codon positions
  # of 'narrow' and 'broad' boundaries
  tb <- split(tb, f = tb$SeqID)
  tb <- lapply(tb, function(ttb){
    boundary_map <- 1:nchar(repertoire[unique(ttb$SeqID)]) -
      stringi::stri_count_fixed(stringi::stri_sub(repertoire[unique(ttb$SeqID)], 1,
                                                  1:nchar(repertoire[unique(ttb$SeqID)])),
                                pattern=".")
    ttb$fiveprime_identical_length <- apply(ttb[, c("start", "fiveprime_identical_length")],
                                            MARGIN = 1, function(x) {
      actual <- boundary_map[x[1]] - x[2]
      o <- which( boundary_map == actual )
      if( length( o ) > 1 ){
        # that must corresponds to a gap. take the outermost pos
        o <- min(o)
      }
      x[1] - o
    })
    ttb$threeprime_identical_length <- apply(ttb[, c("end", "threeprime_identical_length")],
                                             MARGIN = 1, function(x) {
      actual <- boundary_map[x[1]] + x[2]
      o <- which( boundary_map == actual )
      if( length( o ) > 1 ){
        # that must corresponds to a gap. take the outermost pos
        o <- max(o)
      }
      o - x[1]
    })
    ttb
  })
  tb <- do.call("rbind", tb)
  tb$soft_end <- tb$end + tb$threeprime_identical_length
  tb$soft_start <- tb$start - tb$fiveprime_identical_length
  tb
}

#' Give 'narrow' cDNA sequence of gene conversion events in the observed and the germline sequences.
#'
#' @param start numeric, start point of gene conversion event. (column \code{start} from output of \code{batchConvertAnalysis})
#' @param end numeric, end point of gene conversion event. (column \code{start} from output of \code{batchConvertAnalysis})
#' @param start_c numeric, position in the codon (i.e. 1, 2, or 3) for the start point of gene conversion event.
#' @param end_c numeric, position in the codon (i.e. 1, 2, or 3) for the end point of gene conversion event.
#' @param seq_id character, identifier for the sequence of interest. i.e. the relevant entry in the \code{SeqID} column of the output from \code{batchConvertAnalysis}.
#' @param repertoire A named vector of strings storing nucleotide sequences observed from a repertoire. The \code{names} attribute of the vector stores the sequence identifiers.
#' @param seq character, an ungapped DNA sequence of the functional allele.
#'
#' @return A vector with two characters. The first item corresponds to the cDNA sequence corresponding to the location of the annotated gene conversion event, taken from \code{functional}. The second item corresponds to the cDNA sequence at the location of the gene conversion event, observed for the given \code{seq_id} in \code{repertoire}.
#'
getNarrowCDNA <- function( start, end, start_c, end_c,
                           seq_id, repertoire, seq)
{
  start2 <- start - (start_c - 1)
  end2 <- end + (3 - end_c)
  if( any( is.na( c( start, end, start2, end2 ) ) ) ) return(c("", ""))
  germline <- substr( seq, start2, end2 )
  observed <- substr( gsub(".", "", repertoire[seq_id], fixed = TRUE), start2, end2 )
  c(germline, observed)
}

#' Give 'broad' cDNA sequence of gene conversion events in the observed and the germline sequences.
#'
#' @param soft_start numeric, start point of gene conversion event padded with the 5' identical sequence stretch. (column \code{soft_start} from the output of \code{adjustBroadPos}.)
#' @param start numeric, start point of gene conversion event. (column \code{start} from output of \code{adjustBroadPos})
#' @param soft_end numeric, end point of gene conversion event padded with the 3' identical sequence stretch. (column \code{soft_end} from the output of \code{adjustBroadPos}.)
#' @param end numeric, end point of gene conversion event. (column \code{start} from output of \code{adjustBroadPos})
#' @param start_c numeric, position in the codon (i.e. 1, 2, or 3) for the start point of gene conversion event.
#' @param end_c numeric, position in the codon (i.e. 1, 2, or 3) for the end point of gene conversion event.
#' @param seq_id character, identifier for the sequence of interest. i.e. the relevant entry in the \code{SeqID} column of the output from \code{adjustBroadPos}.
#' @param repertoire A named vector of strings storing nucleotide sequences observed from a repertoire. The \code{names} attribute of the vector stores the sequence identifiers.
#' @param seq character, an ungapped DNA sequence of the functional allele.
#'
#' @return A vector with two characters. The first item corresponds to the cDNA sequence corresponding to the 'broad' definition of the annotated gene conversion event, taken from \code{functional}. The second item corresponds to the cDNA sequence for the 'broad' definition of the gene conversion event, observed for the given \code{seq_id} in \code{repertoire}.
#'
getBroadCDNA <- function( soft_start, start, soft_end, end,
                          start_c, end_c, seq_id, repertoire, seq )
{
  start2 <- soft_start - (start_c - 1)
  end2 <- soft_end + (3 - end_c)
  if( any( is.na( c( start, end, start2, end2 ) ) ) ) return(c("", ""))
  germline <- substr( seq, start2, end2 )
  observed <- substr( gsub(".", "", repertoire[seq_id], fixed = TRUE), start2, end2 )
  c(germline, observed)
}

#' Translate amino acid sequences for gene conversion events
#'
#' @param tb data.frame, output from \code{batchConvertAnalysis}.
#' @param repertoire A named vector of strings storing nucleotide sequences observed from a repertoire. The \code{names} attribute of the vector stores the sequence identifiers.
#' @param functional character, filepath to FASTA file containing DNA sequence(s) of the functional V gene allele(s).
#'
#' @return A data.frame with these columns added to \code{tb}:
#' \describe{
#'   \item{germline_AA_narrow}{(Amino acid sequence of the region between \code{start} and \code{end}, from the functional germline allele.}
#'   \item{observed_AA_narrow}{Amino acid sequence of the region between \code{start} and \code{end}, from the observed sequence sampled in the repertoire data.}
#'   \item{germline_AA_broad}{Amino acid sequence of the region between \code{start} and \code{end} plus the flanking stretches given by \code{fiveprime_identical_length} and \code{threeprime_identical_length}, from the functional germline allele.}
#'   \item{observed_AA_broad}{(Amino acid sequence of the region between \code{start} and \code{end} plus the flanking stretches given by \code{fiveprime_identical_length} and \code{threeprime_identical_length}, from the observed sequence sampled in the repertoire data.}
#' }
#'
getAAGeneConversion <- function(tb, repertoire, functional)
{
  tb$threeprime_identical_length <- as.numeric(tb$threeprime_identical_length)
  tb$fiveprime_identical_length <- as.numeric(tb$fiveprime_identical_length)
  # get ungapped coordinates of events for parsing codon positions
  # of 'narrow' and 'broad' boundaries
  tb <- split(tb, f = tb$SeqID)
  tb <- lapply(tb, function(ttb){
    boundary_map <- 1:nchar(repertoire[unique(ttb$SeqID)]) -
      stringi::stri_count_fixed(stringi::stri_sub(repertoire[unique(ttb$SeqID)], 1,
                                                  1:nchar(repertoire[unique(ttb$SeqID)])),
                                pattern=".")
    ttb$start <- sapply(ttb$start, function(x) {
      if(! x %in% 1:nchar(repertoire[unique(ttb$SeqID)])) return(NA)
      boundary_map[x]
    })
    ttb$end <- sapply(ttb$end, function(x) {
      if(! x %in% 1:nchar(repertoire[unique(ttb$SeqID)])) return(NA)
      boundary_map[x]
    })
    ttb
  })
  tb <- do.call("rbind", tb)
  tb$soft_end <- tb$end + tb$threeprime_identical_length
  tb$soft_start <- tb$start - tb$fiveprime_identical_length
  # codon position = remainder(x / 3)
  tb$start_codonpos <- (tb$start %% 3)
  tb$end_codonpos <- (tb$end %% 3)
  tb$start_codonpos <- replace(tb$start_codonpos, which(tb$start_codonpos == 0), 3)
  tb$end_codonpos <- replace(tb$end_codonpos, which(tb$end_codonpos == 0), 3)
  tb$soft_start_codonpos <- (tb$soft_start %% 3)
  tb$soft_end_codonpos <- (tb$soft_end %% 3)
  tb$soft_start_codonpos <- replace(tb$soft_start_codonpos,
                                     which(tb$soft_start_codonpos == 0), 3)
  tb$soft_end_codonpos <- replace(tb$soft_end_codonpos,
                                   which(tb$soft_end_codonpos == 0), 3)
  # get cDNA
  narrow_cdna <- t(apply(tb[, c("start", "end", "start_codonpos", "end_codonpos",
                                 "SeqID", "allele")], MARGIN = 1, function(x){
                                   getNarrowCDNA( as.numeric(x[1]), as.numeric(x[2]), as.numeric(x[3]),
                                                  as.numeric(x[4]), x[5], repertoire = repertoire,
                                                  seq = gsub(".", "", as.character(functional[x[6]]), fixed = TRUE))
                                 }))
  narrow_cdna <- data.frame(narrow_cdna, stringsAsFactors = FALSE)
  colnames(narrow_cdna) <- c("germline_cdna_narrow", "observed_cdna_narrow")
  broad_cdna <- t(apply(tb[, c("soft_start", "start", "soft_end", "end",
                                "soft_start_codonpos", "soft_end_codonpos",
                                "SeqID", "allele")], MARGIN = 1, function(x){
                                  getBroadCDNA( as.numeric(x[1]), as.numeric(x[2]), as.numeric(x[3]),
                                                as.numeric(x[4]), as.numeric(x[5]), as.numeric(x[6]), x[7],
                                                repertoire = repertoire,
                                                seq = gsub(".", "", as.character(functional[x[8]]), fixed = TRUE))
                                }))
  broad_cdna <- data.frame(broad_cdna, stringsAsFactors = FALSE)
  colnames(broad_cdna) <- c("germline_cdna_broad", "observed_cdna_broad")
  # get AA sequence
  narrow_cdna$germline_AA_narrow <-
    suppressWarnings( unname( as.character(Biostrings::translate(
      Biostrings::DNAStringSet( narrow_cdna$germline_cdna_narrow), no.init.codon = TRUE)
    ) ) )
  narrow_cdna$observed_AA_narrow <-
    suppressWarnings( unname( as.character(Biostrings::translate(
      Biostrings::DNAStringSet( narrow_cdna$observed_cdna_narrow), no.init.codon = TRUE)
    ) ) )
  broad_cdna$germline_AA_broad <-
    suppressWarnings( unname( as.character(Biostrings::translate(
      Biostrings::DNAStringSet( broad_cdna$germline_cdna_broad), no.init.codon = TRUE)
    ) ) )
  broad_cdna$observed_AA_broad <-
    suppressWarnings( unname( as.character(Biostrings::translate(
      Biostrings::DNAStringSet( broad_cdna$observed_cdna_broad), no.init.codon = TRUE)
    ) ) )
  broad_cdna[, 3] <- replace(broad_cdna[, 3], which(broad_cdna[, 3] == ""), NA)
  broad_cdna[, 4] <- replace(broad_cdna[, 4], which(broad_cdna[, 4] == ""), NA)
  data.frame( narrow_cdna[, 3:4], broad_cdna[, 3:4], stringsAsFactors = FALSE)
}
