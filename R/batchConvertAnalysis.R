batchConvertAnalysis <- function(functional, pseudogene, repertoire,
                                 blat_exec)
{
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
                       gapwidth = 6)
  toBlat_gw6 <- toBlat_gw6[sapply(toBlat_gw6, function(x) !is.null(x))]
  toBlat_gw6 <- do.call("rbind", toBlat_gw6)
  blat_gw6 <- blat(toBlat_gw6, database = tmpfile_p,
                   blat_exec = blat_exec, min_score = 1)
  blat_gw6 <- formatBlat(blat_gw6, adjustment = 6)

  results_gw6 <- lapply(names(repertoire), ScanGeneConversion, repertoire = repertoire,
                        functional = functional, pseudogenes = pseudogene,
                        blat_all = blat_gw6, blat_whole = blat_whole,
                        lut = functional_mismatches, gapwidth = 6)
  results_gw6 <- results_gw6[sapply(results_gw6, function(x) !is.null(x))]
  results_gw6 <- do.call("rbind", results_gw6)

  # look for unmapped ones and gather a list of their Sequence IDs
  toSecondMap <- unique( results_gw6[which(is.na(results_gw6$gene)), "SeqID"] )

  # second run only for those Sequence IDs with unmapped events. Now merge gaps spaced 3bp apart.
  if( length( toSecondMap ) > 0 ){
    toBlat_gw3 <- lapply(toSecondMap, doBlat,
                         repertoire = repertoire,
                         functional = functional,
                         heavy_or_light = args$heavy_or_light,
                         gapwidth = 3)
    toBlat_gw3 <- toBlat_gw3[sapply(toBlat_gw3, function(x) !is.null(x))]
    toBlat_gw3 <- do.call("rbind", toBlat_gw3)
    blat_gw3 <- blat(toBlat_gw3, database = args$pseudogenes,
                     blat_exec = args$blat, min_score = 1)
    if( is.data.frame(blat_gw3) ){
      if( nrow(blat_gw3) > 0 ){
        blat_gw3 <- formatBlat(blat_gw3, adjustment = 6)
        results_gw3 <- lapply(toSecondMap, ScanGeneConversion, repertoire = repertoire,
                              functional = functional, pseudogenes = pseudogene,
                              blat_all = blat_gw3, lut = functional_mismatches,
                              gapwidth = 3, heavy_or_light = args$heavy_or_light)
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
  results
}
