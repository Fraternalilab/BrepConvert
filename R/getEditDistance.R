#' Function to calculate Levenshtein distance between defined sequence stretches of two longer sequences.
#'
#' @param hits data.frame storing potential hits (donor pseudogenes) that explain a given gene conversion event. Each row corresponds to one definition using a given pseudogene. The boundaries of the pseudogene sequence and the observed repertorie sequence are given.
#' @param query character, the full-length V gene sequence of the observed repertoire sequence for which gene conversion events are annotated.
#' @param subject \code{DNAStringSet} object containing a set of pseudogene DNA sequences.
#' @param method the method used to calculate string distances (default: "levenshtein")
#' 
#' @importFrom Biostrings subseq reverseComplement stringDist
#'
getEditDistance <- function(hits, query, subject, method = "levenshtein")
{
  o <- apply(hits, 1, function(x){
    if (!is.na(x[1])){
      s <- try( Biostrings::subseq(subject[x[1]], as.numeric(x[6]), as.numeric(x[7])), silent = TRUE )
      if(class(s) == "try-error") return(NA)
      if(x[8] == "-"){
        s <- Biostrings::reverseComplement(s)
      }
      query <- gsub(".", "", query, fixed = TRUE)
      q <- Biostrings::subseq(query, as.numeric(x[4]), as.numeric(x[5]))
      as.numeric( Biostrings::stringDist(c(s, q), method = method) )
    } else return(NA)
  })
  names(o) <- hits$subject
  o
}

