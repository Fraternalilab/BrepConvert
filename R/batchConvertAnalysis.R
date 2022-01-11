#' Annotating gene conversion events in BCR repertoire data.
#'
#' @param functional character, filepath to FASTA file containing DNA sequence(s) of the functional V gene allele(s).
#' @param pseudogene character, filepath to FASTA file containing DNA sequence(s) of the pseudogene V gene allele(s).
#' @param repertoire a named vector of characters corresponding to IMGT-gapped DNA sequence from the BCR repertoire. The names are taken as the identifiers of the sequences. See examples below for suggestions on how to generate this from AIRR format repertoire data.
#' @param blat_exec character, filepath to the executable of the BLAT program.
#' @param dist_cutoff vector of 2 integers, 2 distance cut-offs to be used to merge adjacent events. Default = c(6, 3), i.e. the method first try to merge events which are at most 6bp apart to look for pseudogenes donor which can explain this one event; if not, it will try a more stringent cutoff of 3bp.
#' @param convertAA Do you want amino acid sequences of the original, germline functional allele and the observed, converted segment? (default: TRUE)
#'
#' @description This is the main function of the \code{BrepConvert} package for users to annotate gene conversion events in BCR repertoire data, given DNA sequence sets of functional and pseudogene V gene alleles.
#'
#' @return A data.frame with each row corresponding to one gene conversion event. The following annotations are stored in separate columns:
#' \describe{
#'   \item{event}{integer from 1, 2, ... up to n = the number of events observed on a sequence. Just an identifier of gene conversion event.}
#'   \item{possibility}{character from a, b, ... and so on. Denote different possibilities of donor pseudogenes which could account for the observed conversion event.}
#'   \item{start}{integer, position on the repertoire sequence which denotes the start of the conversion event.}
#'   \item{end}{integer, position on the repertoire sequence which denotes the end of the conversion event.}
#'   \item{gene}{character, a semicolon-delimited list of possible donor pseudogenes which could account for the conversion event. \code{NA} if the gene conversion event could not be matched to any pseudogenes.}
#'   \item{fiveprime_identical_length}{integer, the number of nucleotides at the 5' of the named conversion event which is identical between the observed sequence and the named \code{genes}. \code{NA} if the gene conversion event could not be matched to any pseudogenes.}
#'   \item{threeprime_identical_length}{integer, the number of nucleotides at the 3' of the named conversion event which is identical between the observed sequence and the named \code{genes}. \code{NA} if the gene conversion event could not be matched to any pseudogenes.}
#'   \item{edit_distance}{integer, Levenshtein distance comparing the sequence stretch observed on the repertoire sequence and the aligned sequence stretch originated from the donor pseudogene(s). \code{NA} if the gene conversion event could not be matched to any pseudogenes.}
#'   \item{nearest_AID_motif}{integer, the position on the observed sequence where a DNA motif targeted by the AID enzyme can be found closest (at 5') to the gene conversion event.}
#'   \item{AID_motif}{character, the DNA motif targeted by the AID enzyme which is closest (at 5') to the gene conversion event, at the position given by \code{nearest_AID_motif}.}
#'   \item{distance_to_AID_motif}{integer, the number of nucleotides between the named gene conversion event and the given \code{AID_motif}.}
#'   \item{SeqID}{character, identifier for the repertoire sequence, taken from the \code{names} attribute of the input parameter \code{repertoire}.}
#'   \item{seq_event}{character, nucleotide sequence stretch corresponding to the gene conversion event.}
#'   \item{seq_5prime}{character, sequence stretch 10 nucleotides 5' of the start site of the gene conversion event. \code{NA} if the gene conversion event begins at position 1.}
#'   \item{seq_3prime}{character, sequence stretch 10 nucleotides 3' of the end site of the gene conversion event. \code{NA} if the gene conversion event stops at the last position of the V gene.}
#'   \item{germline_AA_narrow}{(avaialble if \code{convertAA == TRUE}) Amino acid sequence of the region between \code{start} and \code{end}, from the functional germline allele.}
#'   \item{observed_AA_narrow}{(avaialble if \code{convertAA == TRUE}) Amino acid sequence of the region between \code{start} and \code{end}, from the observed sequence sampled in the repertoire data.}
#'   \item{germline_AA_broad}{(avaialble if \code{convertAA == TRUE}) Amino acid sequence of the region between \code{start} and \code{end} plus the flanking stretches given by \code{fiveprime_identical_length} and \code{threeprime_identical_length}, from the functional germline allele.}
#'   \item{observed_AA_broad}{(avaialble if \code{convertAA == TRUE}) Amino acid sequence of the region between \code{start} and \code{end} plus the flanking stretches given by \code{fiveprime_identical_length} and \code{threeprime_identical_length}, from the observed sequence sampled in the repertoire data.}
#' }
#'
#' @examples
#' \dontrun{
#' # FASTA files containing pseudogene and functional alleles of the Chicken IGLV
#' # locus are shipped with the package.
#' functional_IGLV <- system.file("extdata/IMGT_Chicken_IGLV_F.fasta",
#'                                package = "BrepConvert")
#' pseudogene_IGLV <- system.file("extdata/IMGT_Chicken_IGLV_P.fasta",
#'                                package = "BrepConvert")
#'
#' # An executable of the BLAT program is also shipped with the package.
#' blat <- system.file("exe/blat", package = "BrepConvert")
#' # NOTE: This works for Linux OS only. For other OS, please download/compile
#' # the executable and indicate the filepath like so:
#' # blat <- "/home/abc/Documents/blat/blat"
#'
#' annotation <- batchConvertAnalysis(
#'   functional = functional_IGLV,
#'   pseudogene = pseudogene_IGLV,
#'   repertoire = repertoire, # notice here the vector is passed, NOT the entire data table!
#'   blat_exec = blat
#' )
#'
#' }
#'
#' @importFrom Biostrings readDNAStringSet writeXStringSet DNAStringSet pairwiseAlignment
#' @importFrom stringr str_extract
#' @importFrom IRanges reduce IRangesList
#'
#' @export
batchConvertAnalysis <- function(functional, pseudogene, repertoire,
                                 blat_exec, dist_cutoff = c(6, 3), convertAA = TRUE)
{
  if( !is.character( repertoire ))
    stop("'repertoire' should be a named vector of characters containing the IMGT-gapped DNA sequences you wish to analyse.")
  if( is.null( names( repertoire ) ) )
    stop("'repertoire' should be a named vector of characters containing the IMGT-gapped DNA sequences you wish to analyse.")
  if( !is.character( functional ) | !is.character( pseudogene ) )
    stop("'functional' and 'pseudogene' should be filepaths pointing to the relevant FASTA files.")
  if( !file.exists( functional ) )
    stop("'functional' should be filepaths pointing to the relevant FASTA files.")
  if( !file.exists( pseudogene ) )
    stop("'pseudogene' should be filepaths pointing to the relevant FASTA files.")
  if( !file.exists( blat_exec ) )
    stop("'blat_exec' should be filepaths pointing to the BLAT executable.")
  blat_msg <- suppressWarnings( system( blat_exec, intern = TRUE ) )
  if( !grepl( "blat - Standalone BLAT", blat_msg[1] ) )
    stop("The BLAT executable does not appear to work. Are you sure you are given run privilege for the executable?")
  if( length( dist_cutoff) != 2 )
    stop("dist_cutoff should be a vector of 2 integers.")
  if( class(dist_cutoff) != "numeric" )
    stop("dist_cutoff should be a vector of 2 integers.")
  dist_cutoff <- as.integer( dist_cutoff )
  dist_cutoff <- sort( dist_cutoff, decreasing = TRUE )
  if( !is.logical( convertAA ))
    stop("'convertAA' must be either TRUE or FALSE.")

  seqnames <- names(repertoire)
  # read in functional gene
  functional <- Biostrings::readDNAStringSet(functional)
  pseudogene <- Biostrings::readDNAStringSet(pseudogene)
  names(functional) <- stringr::str_extract(names(functional), "IG.V[0-9A-Z\\-\\*]*\\||IG.V.*$")
  names(pseudogene) <- stringr::str_extract(names(pseudogene), "IG.V[0-9A-Z\\-\\*]*\\||IG.V.*$")
  names(functional) <- gsub("|", "", names(functional), fixed = TRUE)
  names(pseudogene) <- gsub("|", "", names(pseudogene), fixed = TRUE)

  # write out an ungapped version of pseudogenes to tempfile for blat
  tmpfile_p <- tempfile()
  Biostrings::writeXStringSet(
    Biostrings::DNAStringSet(sapply(pseudogene, function(x) gsub(".", "", x, fixed = TRUE))),
    filepath = tmpfile_p
  )
  pseudogene <- Biostrings::readDNAStringSet(tmpfile_p)

  # Step (0) Generate look-up tables (IRanges)
  # 0a. Look-up table of AID hotspot locations
  #functional_AID <- sort(Reduce(c, sapply(c("AGC", "AGT", "GGC", "GGT", "TGC"),
  #                                        function(x) vmatchPattern(x, subject = functional)[[1]])))

  # 0b. Look-up table of mismatches of each pseudogene compared against the functional allele
  functional_mismatches <- lapply(functional, function(f){
    m <- Biostrings::pairwiseAlignment(subject = f, pseudogene, type = "local")
    m <- sapply(m, function(x){
      m <- IRanges::reduce(c(as(x@pattern@mismatch, "IRangesList")[[1]], x@pattern@indel[[1]]))
      if(length(m) > 0) names(m) <- rep(names(x@pattern@unaligned), length(m))
      m
    })
    m <- Reduce( c, m )
  })
  names(functional_mismatches) <- names(functional)

  # blat of whole sequence; to use if BLATing the short stretch doesn't yield results
  blat_whole <- blat(data.frame(qname = seqnames,
                                seq = gsub(".", "", repertoire, fixed = TRUE),
                                stringsAsFactors = FALSE),
                     database = tmpfile_p, blat_exec = blat_exec)
  blat_whole[, "gene"] <- gsub("|", "",
                               stringr::str_extract(blat_whole[, "gene"],
                                                    "IG.V[0-9A-Z\\-\\*]*\\||IG.V.*$"),
                               fixed = TRUE)
  blat_whole <- blat_whole[, c("qname", "gene", "mismatch", "gaps", "q_start", "q_end",
                               "s_start", "s_end")]
  if( is.data.frame(blat_whole) ){
    if( nrow(blat_whole) > 0){
      colnames(blat_whole) <- c("query", "subject", "mismatches", "gaps_open", "q_start",
                                "q_end", "s_start", "s_end")
    } else blat_whole <- data.frame()
  } else blat_whole <- data.frame()

  # first run: merge gaps spaced 6bp apart
  toBlat_gw6 <- lapply(names(repertoire), doBlat,
                       repertoire = repertoire,
                       functional = functional,
                       gapwidth = dist_cutoff[1])
  toBlat_gw6 <- toBlat_gw6[sapply(toBlat_gw6, function(x) !is.null(x))]
  toBlat_gw6 <- do.call("rbind", toBlat_gw6)
  blat_gw6 <- blat(toBlat_gw6, database = tmpfile_p,
                   blat_exec = blat_exec, min_score = 1)
  blat_gw6 <- formatBlat(blat_gw6, adjustment = dist_cutoff[1])

  results_gw6 <- lapply(names(repertoire), ScanGeneConversion, repertoire = repertoire,
                        functional = functional, pseudogenes = pseudogene,
                        blat_all = blat_gw6, blat_whole = blat_whole,
                        lut = functional_mismatches, gapwidth = dist_cutoff[1])
  results_gw6 <- results_gw6[sapply(results_gw6, function(x) !is.null(x))]
  results_gw6 <- do.call("rbind", results_gw6)

  # look for unmapped ones and gather a list of their Sequence IDs
  toSecondMap <- unique( results_gw6[which(is.na(results_gw6$gene)), "SeqID"] )

  # second run only for those Sequence IDs with unmapped events. Now merge gaps spaced 3bp apart.
  if( length( toSecondMap ) > 0 ){
    toBlat_gw3 <- lapply(toSecondMap, doBlat,
                         repertoire = repertoire,
                         functional = functional,
                         gapwidth = dist_cutoff[2])
    toBlat_gw3 <- toBlat_gw3[sapply(toBlat_gw3, function(x) !is.null(x))]
    toBlat_gw3 <- do.call("rbind", toBlat_gw3)
    blat_gw3 <- try( blat(toBlat_gw3, database = tmpfile_p,
                          blat_exec = blat_exec, min_score = 1), silent = TRUE )
    if( is.data.frame(blat_gw3) ){
      if( nrow(blat_gw3) > 0 ){
        blat_gw3 <- formatBlat(blat_gw3, adjustment = dist_cutoff[2])
        results_gw3 <- lapply(toSecondMap, ScanGeneConversion, repertoire = repertoire,
                              functional = functional, pseudogenes = pseudogene,
                              blat_all = blat_gw3, blat_whole = blat_whole,
			      lut = functional_mismatches, gapwidth = dist_cutoff[2])
        results_gw3 <- results_gw3[sapply(results_gw3, function(x) !is.null(x))]
        results_gw3 <- do.call("rbind", results_gw3)
      } else results_gw3 <- data.frame()
    } else results_gw3 <- data.frame()
  } else {
    results_gw3 <- data.frame()
  }

  if( nrow(results_gw3) > 0 ){
    # merge the two runs together with 6bp run having precedence.
    results <- lapply(unique(results_gw6$SeqID), function(x){
      gw3 <- results_gw3[ which(results_gw3$SeqID == x), ]
      gw6 <- results_gw6[ which(results_gw6$SeqID == x), ]
      if( nrow(gw3) > 0 & nrow(gw6) > 0){
        # look for (smaller) matches in the gw3 scan and replace the hits in gw6
        # with gene NA (ie events without matched pseudogenes)
        to_add <- c()
        to_delete <- c()
        for(i in 1:nrow(gw6)){
          if( is.na(gw6[i, "gene"]) ){
            a <- which(sapply(1:nrow(gw3), function(j){
              intersection <- intersect( IRanges::IRanges(gw6[i, "start"], gw6[i, "end"]),
                                         IRanges::IRanges(gw3[j, "start"], gw3[j, "end"]))
              length(intersection) > 0
            }))
            if(length(a) > 0){
              to_add <- c(to_add, a)
              to_delete <- c(to_delete, i)
            }
          }
        }
        if(length(to_delete) > 0) gw6 <- gw6[-to_delete, ]
        if(length(to_add) > 0) gw6 <- rbind(gw6, gw3[to_add, ])
        # gw6 <- gw6[order(gw6$start, gw6$edit_distance, decreasing = FALSE), ]
      }
      # number events and possibility
      gw6 <- split(gw6, f = gw6[, c("start")], drop = TRUE)
      gw6 <- do.call("rbind", lapply(1:length(gw6), function(k){
        tb <- gw6[[k]]
        tb$event <- k
        tb$possibility <- letters[1:nrow(tb)]
        tb[, c( ncol(tb)-1, ncol(tb), 1:(ncol(tb) - 2) )]
      }))
      return(gw6)
    })
  } else {
    results <- lapply(unique(results_gw6$SeqID), function(x){
      gw6 <- results_gw6[ which(results_gw6$SeqID == x), ]
      # number events and possibility
      gw6 <- split(gw6, f = gw6[, c("start")], drop = TRUE)
      gw6 <- do.call("rbind", lapply(1:length(gw6), function(k){
        tb <- gw6[[k]]
        tb$event <- k
        tb$possibility <- letters[1:nrow(tb)]
        tb[, c( ncol(tb)-1, ncol(tb), 1:(ncol(tb) - 2) )]
      }))
      return(gw6)
    })
  }
  results <- do.call("rbind", results)
  if( convertAA ){
    results <- cbind( results, getAAGeneConversion( results, repertoire,
                                                    functional = functional) )
  }
  results[, -which(colnames(results) == "allele")]
}
