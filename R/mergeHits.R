mergeHits <- function(hits, lut_functional)
{
  # merge consecutive hits on the same pseudogene
  n = 1
  hits$event <- n
  genes <- lapply(hits$gene, function(z) {
    if(!is.na(z)) unlist(strsplit(z, split = ";")) else z
  })
  # only one event
  if(nrow(hits) == 1) {
    if( is.na(hits[1, "gene"]) ){
      return(data.frame(#event = 1, possibility = "a",
        start = hits[1, "start"],
        end = hits[1, "end"], gene = NA,
        fiveprime_identical_length = NA,
        threeprime_identical_length = NA,
        edit_distance = NA, stringsAsFactors = FALSE
      ))
    } else {
      edits <- unlist(strsplit(hits[1, "edit_distance"], split = ";"))
      fp <- unlist(strsplit(hits[1, "fiveprime_identical_length"], split = ";"))
      tp <- unlist(strsplit(hits[1, "threeprime_identical_length"], split = ";"))
      names(edits) <- genes[[1]]
      names(fp) <- genes[[1]]
      names(tp) <- genes[[1]]
      edits <- edits[ which(edits == min(edits, na.rm = TRUE)) ]
      fp <- fp[ names(edits) ]
      tp <- tp[ names(edits) ]
      if(length(fp) == 0) {
        hits <- data.frame(start = hits[1, "start"], end = hits[1, "end"],
                           gene = NA, fiveprime_identical_length = NA,
                           threeprime_identical_length = NA,
                           edit_distance = NA, event = 1,
                           start_ungapped = hits[1, "start_ungapped"],
                           end_ungapped = hits[1, "end_ungapped"],
                           stringsAsFactors = FALSE)
      } else {
        hits <- data.frame(start = hits[1, "start"], end = hits[1, "end"],
                           gene = names(fp), fiveprime_identical_length = fp,
                           threeprime_identical_length = tp,
                           edit_distance = edits, event = 1,
                           start_ungapped = hits[1, "start_ungapped"],
                           end_ungapped = hits[1, "end_ungapped"],
                           stringsAsFactors = FALSE)
      }
      hits <- ddply(hits, .variables = c("event", "fiveprime_identical_length",
                                         "threeprime_identical_length" ,"edit_distance"),
                    summarise,
                    start = min(start), end = max(end),
                    start_ungapped = min(start_ungapped),
                    end_ungapped = max(end_ungapped),
                    gene = paste(gene, collapse = ";")
      )
      return( hits[, c("start", "end", "gene",
                       "fiveprime_identical_length", "threeprime_identical_length",
                       "edit_distance", "start_ungapped", "end_ungapped")] )
    }
  }
  # more than one event
  for(z in 2:nrow(hits)){
    matches <- match(genes[[z]],
                     try( unlist(strsplit(hits[z-1, "gene"], split = ";")),
                          silent = TRUE) )
    if(length(which(!is.na(matches))) > 0){
      # count mismatches (against germline) between hits; if 0, merge; otherwise leave it
      o <- sapply(genes[[z]][ which(!is.na(matches)) ], function(g){
        spacer_start <- as.numeric( unlist(strsplit(hits[z-1, "q_end"], split = ";")) ) + 1
        names(spacer_start) <- genes[[z-1]]
        spacer_end <- as.numeric( unlist(strsplit(hits[z, "q_start"], split = ";")) ) - 1
        names(spacer_end) <- genes[[z]]
        stretch <- try( IRanges::IRanges(start = spacer_start[g], end = spacer_end[g]), silent = TRUE )
        if( class(stretch) == "try-error" ) return(1000) # just an arbitrary large number
        else{
          m <- sum(
            IRanges::width(
              IRanges::intersect(stretch,
                                 lut_functional[ which(names(lut_functional) == g) ])
            ))
          m
        }
      })
      names(o) <- genes[[z]][ which(!is.na(matches)) ]
      o <- o[o == 0]
      if(length(o) > 0) {
        # sum over edit distances
        e_z1 <- as.numeric( unlist(strsplit(hits[z-1, "edit_distance"], split = ";")) )
        names(e_z1) <- unlist(strsplit(hits[z-1, "gene"], split = ";"))
        e_z <- as.numeric( unlist(strsplit(hits[z, "edit_distance"], split = ";")) )
        names(e_z) <- genes[[z]]
        e_z1 <- e_z1[names(o)]; e_z <- e_z[names(o)]
        edits <- sapply(names(o), function(g) e_z1[g] + e_z[g])
        names(edits) <- names(o)
        # select the gene with min edit distance
        edits <- edits[ which(edits == min(edits, na.rm = TRUE)) ]
        hits[z, "edit_distance"] <- paste(edits, collapse = ";")
        hits[z-1, "edit_distance"] <- paste(edits, collapse = ";")
        # copy over fiveprime identity stretch length
        fp <- unlist(strsplit(hits[z-1, "fiveprime_identical_length"], split = ";"))
        names(fp) <- unlist(strsplit(hits[z-1, "gene"], split = ";"))
        hits[z, "fiveprime_identical_length"] <- paste(fp[names(edits)], collapse = ";")
        hits[z-1, "fiveprime_identical_length"] <- paste(fp[names(edits)], collapse = ";")
        o <- paste(names(edits), collapse = ";")
        hits[z, "gene"] <- o; hits[z-1, "gene"] <- o
        hits[z, "event"] <- hits[z-1, "event"]
        hits[z, "start"] <- hits[z-1, "start"]
        hits[z, "start_ungapped"] <- hits[z-1, "start_ungapped"]
        hits[z-1, "end"] <- hits[z, "end"]
        hits[z-1, "end_ungapped"] <- hits[z, "end_ungapped"]
      } else {
        n <- n + 1
        hits[z, "event"] <- n
      }
    }
    else {
      if( hits[z, "start"] != hits[z-1, "start"] ){
        n <- n + 1
      }
      hits[z, "event"] <- n
    }
  }
  hits <- hits[, c("start", "end", "gene", "fiveprime_identical_length",
                   "threeprime_identical_length" , "edit_distance", "event",
                   "start_ungapped", "end_ungapped")]
  hits <- unique(hits)
  to_delete <- c()
  if(nrow(hits) >= 2){
    for(i in 2:nrow(hits)){
      if( hits[i, "start"] == hits[i-1, "start"] & hits[i, "end"] != hits[i-1, "end"] ){
        to_delete <- c(to_delete, i-1)
      }
    }
  }
  if(length(to_delete) > 0) hits <- hits[-to_delete, ]
  # choose gene with min edit distance
  g <- apply(hits[, c("gene", "edit_distance", "fiveprime_identical_length", "threeprime_identical_length")], MARGIN = 1, function(x){
    if(is.na(x[1]) | is.na(x[2])) return( list("edit" = NA, "fp" = NA, "tp" = NA ) )
    else{
      genes <- unlist(strsplit(x[1], split = ";"))
      edits <- as.numeric( unlist(strsplit(x[2], split = ";")) )
      names(edits) <- genes
      fp <- unlist(strsplit(x[3], split = ";"))
      names(fp) <- genes
      edits <- edits[ which(edits == min(edits, na.rm = TRUE)) ]
      fp <- fp[ names(edits) ]
      tp <- unlist(strsplit(x[4], split = ";"))
      names(tp) <- genes
      edits <- edits[ which(edits == min(edits, na.rm = TRUE)) ]
      tp <- tp[ names(edits) ]
      list("edit" = edits, "fp" = fp, "tp" = tp)
    }
  })
  hits$gene <- sapply(g, function(h) paste(names(h$edit), collapse = ";"))
  hits$edit_distance <- sapply(g, function(h) paste(h$edit, collapse = ";"))
  hits$fiveprime_identical_length <- sapply(g, function(h) paste(h$fp, collapse = ";"))
  hits$threeprime_identical_length <- sapply(g, function(h) paste(h$tp, collapse = ";"))
  # deconvolute ";"-separated into separated rows
  hits <- split(hits, hits$event)
  hits <- do.call("rbind", lapply(hits, function(tb){
    if( is.na(tb[1, "gene"]) ){
      o <- data.frame(start = tb[1, "start"], end = tb[1, "end"],
                      gene = NA, fiveprime_identical_length = NA,
                      threeprime_identical_length = NA,
                      edit_distance = NA, event = tb[1, "event"],
                      start_ungapped = tb[1, "start_ungapped"],
                      end_ungapped = tb[1, "end_ungapped"],
                      stringsAsFactors = FALSE)
    } else if( tb[1, "gene"] == "" ){
      o <- data.frame(start = tb[1, "start"], end = tb[1, "end"],
                      gene = NA, fiveprime_identical_length = NA,
                      threeprime_identical_length = NA,
                      edit_distance = NA, event = tb[1, "event"],
                      start_ungapped = tb[1, "start_ungapped"],
                      end_ungapped = tb[1, "end_ungapped"],
                      stringsAsFactors = FALSE)
    } else {
      o <- data.frame(start = tb[1, "start"], end = tb[1, "end"],
                      gene = unlist(strsplit(tb[1, "gene"], split = ";")),
                      fiveprime_identical_length = unlist(strsplit(tb[1, "fiveprime_identical_length"], split = ";")),
                      threeprime_identical_length = unlist(strsplit(tb[1, "threeprime_identical_length"], split = ";")),
                      edit_distance = unlist(strsplit(tb[1, "edit_distance"], split = ";")),
                      event = tb[1, "event"],
                      start_ungapped = tb[1, "start_ungapped"],
                      end_ungapped = tb[1, "end_ungapped"],
                      stringsAsFactors = FALSE)
    }
    o
  }))
  hits <- ddply(hits, .variables = c("event", "fiveprime_identical_length", "threeprime_identical_length", "edit_distance"),
                summarise,
                start = min(start), end = max(end),
                start_ungapped = min(start_ungapped), end_ungapped = max(end_ungapped),
                gene = paste(gene, collapse = ";")
  )
  hits[, c("start", "end", "gene",
           "fiveprime_identical_length", "threeprime_identical_length",
           "edit_distance", "start_ungapped", "end_ungapped")]
}

