getIdenticalLength <- function(hits, subject, query, type, max_dist = 100)
{
  # test (up to 'max') whether the 5' stretches of the subject and the query are identical
  o <- apply(hits, 1, function(x){
    oo <- sapply(1:max_dist, function(w){
      if(w >= as.numeric(x[2])) return(FALSE)
      if( (x[8] == "+" & type == '5p') ){
        if(w >= as.numeric(x[6])) return(FALSE)
        s <- try( subseq(subject[x[1]], as.numeric(x[6]) - w, as.numeric(x[6])), silent = TRUE )
      } else if( (x[8] == "-" & type == '3p') ){
        if(w >= as.numeric(x[4])) return(FALSE)
        s <- try( reverseComplement(subseq(subject[x[1]], as.numeric(x[6]) - w, as.numeric(x[6]))), silent = TRUE )
      } else if((x[8] == "-" & type == '5p') ){
        if(w >= as.numeric(x[5])) return(FALSE)
        s <- try( reverseComplement(subseq(subject[x[1]], as.numeric(x[7]), as.numeric(x[7]) + w)),
                  silent = TRUE )
      } else if( (x[8] == "+" & type == '3p') ){
        if(w >= as.numeric(x[5])) return(FALSE)
        s <- try( subseq(subject[x[1]], as.numeric(x[7]), as.numeric(x[7]) + w), silent = TRUE )
      }
      if(class(s) == "try-error") return(NA)
      query <- gsub(".", "", query, fixed = TRUE)
      if (type == '5p') {
        q <- try( subseq( query, as.numeric(x[4]) - w, as.numeric(x[4]) ), silent = TRUE )
      } else if(type == '3p') {
        q <- try( subseq( query, as.numeric(x[5]), as.numeric(x[5]) + w ), silent = TRUE )
      }
      if( class(q) == "try-error" ) return(NA)
      return(s == q)
    })
    oo <- oo[!is.na(oo)]
    if(length(oo) == 0) return(NA) else return( sum(oo) )
  })
  names(o) <- hits$subject
  o
}

