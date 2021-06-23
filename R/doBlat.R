doBlat <- function(seqname, repertoire, functional,
                   gapwidth = 3, flank=6)
{
  # Step (1) Alignment and identify stretches to be mapped
  # 1a. alignment
  # aln to every allele in the functional list. Take the allele with higher Sequence Identity.
  aln <- lapply(functional, function(f){
    F_alignment <- pairwiseAlignment(f, repertoire[seqname], type="global-local", gapOpening=20)#type = "local")
    F_alignment
  })
  allele <- names(aln)[which.max(sapply(aln, pid))]
  F_alignment <- aln[[which.max(sapply(aln, pid))]]
  # 1b. if Percent ID is below 50 - likely it is an incomplete read; skip this.
  # If not aligned from the beginning of the leader, to be safe, filter away
  if( pid(F_alignment) < 50 ) return(NULL)
  if( start(F_alignment@pattern@range) != 1 ) return(NULL)

  # 1c. process mismatch regions
  mismatches <- as(F_alignment@pattern@mismatch, "IRangesList")
  mismatches <- mismatches[[1]]
  # group together mismatches if they are less than 6bp apart
  mismatches <- reduce(mismatches, min.gapwidth=gapwidth)
  # stretches at least 3bp are 'gaps' to be matched to pseudogenes
  gaps <- mismatches[width(mismatches) >= 3]
  # gaps_flanked <- reduce(gaps + 6) # add flanking stretch to help alignment
  # those less than 3bp are mismatches occur outside of gene-converted regions; to be reported
  mismatches <- mismatches[width(mismatches) < 3]

  # if no gaps, skip
  if(length(gaps) == 0) return(NULL)

  # remove the IMGT gaps for now to fetch flanking bits
  ungapped_sequence <- gsub(".", "", repertoire[seqname], fixed = TRUE)

  # if aln doesn't begin at position 1 pad the gap ranges accordingly
  if( start( F_alignment@pattern@range ) != 1 ){
    gaps <- shift( gaps, start( F_alignment@pattern@range ) - 1 )
  }

  # if beginning of repertoire sequence is gapped, ignore the gap begins at pos 1
  if( grepl("^\\.", repertoire[seqname]) ){
    gaps <- gaps[start(gaps) > 1]
  }
  # adjust the positions to give start and end points in the ungapped version
  gaps <- as.data.frame( gaps )
  gaps$start_ungapped <- sapply(gaps[, 1], function(x){
    n_dots <- stringr::str_count(substr(repertoire[seqname], start = 1, stop = x),
                                 pattern = stringr::fixed("."))
    x - n_dots
  })
  gaps$end_ungapped <- sapply(gaps[, 2], function(x){
    n_dots <- stringr::str_count(substr(repertoire[seqname], start = 1, stop = x),
                                 pattern = stringr::fixed("."))
    x - n_dots
  })
  do.call("rbind", apply(gaps, 1, function(x){
    data.frame(qname = paste0(seqname, "/", x[1], "-", x[2]),
               seq = as.character(subseq(ungapped_sequence, max(1, x[4]-flank),
                                         min(nchar(ungapped_sequence), x[5]+flank))),
               stringsAsFactors = FALSE)
  }))
}
