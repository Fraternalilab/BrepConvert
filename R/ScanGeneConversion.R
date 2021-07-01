ScanGeneConversion <- function(seqname, repertoire, functional,
                               pseudogenes, blat_all, blat_whole, lut,
                               gapwidth = 3){
  cat(paste0(seqname, ' ...\n'))
  # Step (1) Alignment and identify stretches to be mapped
  # 1a. alignment
  # aln to every allele in the functional list. Take the allele with higher Sequence Identity.
  aln <- lapply(functional, function(f){
    F_alignment <- Biostrings::pairwiseAlignment(f, repertoire[seqname], type="global-local", gapOpening=20)#type = "local")
    F_alignment
  })
  allele <- names(aln)[which.max(sapply(aln, Biostrings::pid))]
  F_alignment <- aln[[which.max(sapply(aln, Biostrings::pid))]]
  # 1b. if Percent ID is below 50 - likely it is an incomplete read; skip this.
  # If not aligned from the beginning of the leader, to be safe, filter away
  if( Biostrings::pid(F_alignment) < 50 ) return(NULL)
  if( IRanges::start(F_alignment@pattern@range) != 1 ) return(NULL)

  # 1c. process mismatch regions
  mismatches <- as(F_alignment@pattern@mismatch, "IRangesList")
  mismatches <- mismatches[[1]]
  # group together mismatches if they are less than 6bp apart
  mismatches <- IRanges::reduce(mismatches, min.gapwidth=gapwidth)
  # stretches at least 3bp are 'gaps' to be matched to pseudogenes
  gaps <- mismatches[IRanges::width(mismatches) >= 3]
  # gaps_flanked <- reduce(gaps + 6) # add flanking stretch to help alignment
  # those less than 3bp are mismatches occur outside of gene-converted regions; to be reported
  mismatches <- mismatches[IRanges::width(mismatches) < 3]
  # if beginning of repertoire sequence is gapped, ignore the gap begins at pos 1
  if( grepl("^\\.", repertoire[seqname]) ){
    gaps <- gaps[IRanges::start(gaps) > 1]
  }
  # if aln doesn't begin at position 1 pad the gap ranges accordingly
  gaps <- IRanges::shift( gaps, IRanges::start( F_alignment@pattern@range ) - 1 )

  # if no gaps, skip
  if(length(gaps) == 0) return(NULL)

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


  # Step (2) Identify pseudogenes
  # 2a. process blat; identify top 5 hits; skip if no match
  # (i) blat with substrings
  if( is.data.frame( blat_all ) ){
    if( nrow(blat_all) > 0 ){
      blat_part <- blat_all[blat_all[, 1] == seqname, ]
    } else blat_part <- data.frame()
  } else blat_part <- data.frame()
  if(nrow(blat_part) == 0) return(NULL)

  # (ii) blat with whole sequence
  if( is.data.frame( blat_whole ) ){
    if( nrow(blat_whole) > 0 ){
      blat_whole <- blat_whole[blat_whole[, 1] == seqname, ]
    } else blat_whole <- data.frame()
  } else blat_whole <- data.frame()
  if(nrow(blat_whole) == 0) return(NULL)


  # 2b. alignment only of the blat-matched regions. Check edit distance of the gaps
  # Report top pseudogene match and edit distance.

  # this is the amount of shift to convert from sequence starting at leader to
  # numbering starting from the coding segment

  top_hits <- apply(gaps, MARGIN = 1, function(x){
    hit_range <- IRanges::IRanges(start = x[1], end = x[2])
    startpos <- IRanges::start(hit_range); endpos <- IRanges::end(hit_range)
    matches <- getTopBlat(blat_part, start = startpos, end = endpos, ntop = 1,
                          adjustment = NULL)
    if(nrow(matches) == 0){
      # try looking at the blat result on whole sequence
      # note the boundaries weren't corrected here so we use the ungapped positions
      matches <- getTopBlat(blat_whole, start = x[4], end = x[5], ntop = NULL,
                            adjustment = NULL)
      startpos <- x[4]; endpos <- x[5]
      if(nrow(matches) == 0) {
        # then call this off and go to the next
        return( data.frame(start = startpos, end = endpos,
                           q_start = NA, q_end = NA, gene = NA,
                           fiveprime_identical_length = NA,
                           threeprime_identical_length = NA,
                           edit_distance = NA, #direction = NA,
                           start_ungapped = x[4], end_ungapped = x[5],
                           stringsAsFactors = FALSE) )
      }
    }
    matches <- do.call("rbind", apply(matches, 1, function(y){
      if(as.numeric(y[3]) > as.numeric(y[2])) dir = '+' else dir = '-'
      if(as.numeric(y[3]) < as.numeric(y[2])){
        # in reverse
        s_start <- as.numeric(y[3]) - abs(startpos - as.numeric(y[4]))
        s_end <- s_start - (x[5] - x[4])
        data.frame(subject = y[1], q_start = startpos, q_end = endpos,
                   q_start_ungapped = x[4], q_end_ungapped = x[5],
                   s_start = s_start, s_end = s_end,
                   dir = dir, stringsAsFactors = FALSE)
      } else {
        s_start <- as.numeric(y[2]) + abs(startpos - as.numeric(y[4]))
        s_end <- s_start + (x[5] - x[4])
        data.frame(subject = y[1], q_start = startpos, q_end = endpos,
                   q_start_ungapped = x[4], q_end_ungapped = x[5],
                   s_start = s_start, s_end = s_end,
                   dir = dir, stringsAsFactors = FALSE)
      }
    }))
    edits <- getEditDistance(matches, query = repertoire[seqname],
                             subject = pseudogenes)
    fiveprime <- getIdenticalLength(matches, query = repertoire[seqname],
                                    subject = pseudogenes, type = '5p')
    threeprime <- getIdenticalLength(matches, query = repertoire[seqname],
                                     subject = pseudogenes, type = '3p')
    #    edits <- as.data.frame(edits); rownames(edits) <- matches$subject
    #    fiveprime <- as.data.frame(fiveprime); rownames(fiveprime) <- matches$subject
    d <- merge(as.data.frame(fiveprime), as.data.frame(threeprime), by = "row.names")
    colnames(d)[1] <- 'gene'
    d <- merge(d, as.data.frame(edits), by.x = 'gene', by.y = "row.names")
    # d <- as.data.frame(fiveprime)
    matches <- merge(matches, d, by.x = "subject", by.y = "gene")
    #matches <- ddply(matches, .variables = c("fiveprime", #"edits",
    #                                         "dir"), summarise,
    #                 gene = paste(subject, collapse = ";"))

    # if gap width > 10, reject solutions where edit_distance >= 0.5 * width
    if(IRanges::width(hit_range) > 10){
      matches <- matches[which(matches$edits < 0.5 * (endpos - startpos + 1)), ]
    }
    matches <- matches[order(matches$fiveprime, matches$threeprime, decreasing = TRUE), ]
    # matches$fiveprime <- replace(matches$fiveprime, which(matches$fiveprime == 20), ">=20")
    colnames(matches)[1] <- 'gene'
    unique(data.frame(start = x[1], end = x[2],
                      q_start = paste(matches$s_start, collapse = ";"),
                      q_end = paste(matches$s_end, collapse = ";"),
                      gene = paste(matches$gene, collapse = ";"),
                      fiveprime_identical_length = paste(matches$fiveprime, collapse = ";"),
                      threeprime_identical_length = paste(matches$threeprime, collapse = ";"),
                      edit_distance = paste(matches$edits, collapse = ";"),
                      start_ungapped = x[4], end_ungapped = x[5],
                      stringsAsFactors = FALSE))
  })
  top_hits <- top_hits[sapply(top_hits, function(x) !is.null(x))]
  if(is.null(top_hits)) return(NULL)
  top_hits <- do.call("rbind", top_hits)
  # Merge consecutive events involving the same pseudogene
  top_hits <- mergeHits(top_hits, lut_functional = lut[[allele]])
  top_hits$SeqID <- seqname
  # extract substrings
  ungapped_sequence <- gsub(".", "", repertoire[seqname], fixed = TRUE)
  top_hits$seq_event <- apply(top_hits[, c("start_ungapped", "end_ungapped")],
                              MARGIN = 1, function(x){
    toupper( as.character( Biostrings::subseq( ungapped_sequence, x[1], x[2]) ) )
  })

  # get nearest AID motif
  AID <- Reduce(c,
                sapply(c("AGC", "AGT", "GGC", "GGT", "TGC", "AAC", "TAC"), function(x) {
                  o <- Biostrings::vmatchPattern(pattern = x, subject = functional[allele])[[1]]
                  S4Vectors::values(o) <- S4Vectors::DataFrame(motif = x)
                return(o)
         })
  )#repertoire[seqname])[[1]])))
  #AID <- AID[start(AID) > 58]
  AID <- AID[IRanges::width(AID) >= 3]
  #AID <- shift(AID, -58 + 1)

  # filter away identified gaps and mismatches
  #AID <- setdiff(setdiff(AID, mismatches), IRanges(start = top_hits$start, end = top_hits$end))
  # AID <- AID[(AID %outside% mismatches)]
  # AID <- AID[(AID %outside% IRanges(start = top_hits$start, end = top_hits$end))]

  top_hits$nearest_AID_motif <- apply(top_hits, MARGIN = 1, function(x){
    # genes <- unlist(strsplit(x[5], split = ";"))
    gap <- IRanges::IRanges(as.numeric(x[1]), as.numeric(x[2]))
    nearest_index <- IRanges::follow(gap, AID)
    if(!is.na(nearest_index)){
      IRanges::start(AID[nearest_index])
    } else {
      # then the first AID motif must coincide the 'start' of the gap
      IRanges::start(AID[1])
    }
  })

  top_hits$AID_motif <- apply(top_hits, MARGIN = 1, function(x){
    # genes <- unlist(strsplit(x[5], split = ";"))
    gap <- IRanges::IRanges(as.numeric(x[1]), as.numeric(x[2]))
    nearest_index <- IRanges::follow(gap, AID)
    if(!is.na(nearest_index)){
      S4Vectors::values(AID[nearest_index])[1, 1]
    } else {
      # then the first AID motif must coincide the 'start' of the gap
      S4Vectors::values(AID[1])[1, 1]
    }
  })

  top_hits$distance_to_AID_motif <- apply(top_hits, MARGIN = 1, function(x){
    # genes <- unlist(strsplit(x[5], split = ";"))
    gap <- IRanges::IRanges(as.numeric(x[1]), as.numeric(x[2]))
    nearest_index <- IRanges::follow(gap, AID)
    if(!is.na(nearest_index)){
      IRanges::distance(gap, AID[nearest_index])
    } else {
      # then the first AID motif must coincide the 'start' of the gap
      0
    }
    # distance(gap, AID[follow(gap, AID)])
  })

  # extract 5' and 3' strings relative to the conversion events
  top_hits$seq_5prime <- sapply(top_hits[, "start"], function(x){
    if( (max(1, x - 10) < (x - 1)) & (max(1, x - 10) < nchar(as.character(functional[allele]))) ){
      s <- as.character( Biostrings::subseq( functional[allele], 1,
                                             min(nchar(as.character(functional[allele])), x - 1)) )
      s <- gsub(".", "", s, fixed = TRUE)
      Biostrings::subseq(s, max(1, nchar(s) - 9), nchar(s))
    } else return(NA)
  })
  top_hits$seq_3prime <- sapply(top_hits[, "end"], function(x){
    if( ((x + 1) <  min(nchar(as.character(functional[allele])), x + 10)) &
        (x + 1 < nchar(as.character(functional[allele]))) ){
      s <- as.character( Biostrings::subseq( functional[allele], max(1, x + 1),
                                             nchar(as.character(functional[allele]))) )
      s <- gsub(".", "", s, fixed = TRUE)
      Biostrings::subseq(s, 1, min(10, nchar(s)))
    } else return(NA)
  })

  top_hits[, c("start","end","gene","fiveprime_identical_length",
               "threeprime_identical_length","edit_distance",
               "nearest_AID_motif","AID_motif","distance_to_AID_motif",
               "SeqID","seq_event","seq_5prime","seq_3prime")]
}
