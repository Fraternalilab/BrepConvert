getEditDistance <- function(hits, query, subject, method = "levenshtein")
{
  o <- apply(hits, 1, function(x){
    if (!is.na(x[1])){
      s <- try( subseq(subject[x[1]], as.numeric(x[6]), as.numeric(x[7])), silent = TRUE )
      if(class(s) == "try-error") return(NA)
      if(x[8] == "-"){
        s <- reverseComplement(s)
      }
      query <- gsub(".", "", query, fixed = TRUE)
      q <- subseq(query, as.numeric(x[4]), as.numeric(x[5]))
      as.numeric( Biostrings::stringDist(c(s, q), method = method) )
    } else return(NA)
  })
  names(o) <- hits$subject
  o
}

