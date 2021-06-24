getIdenticalLength <- function(hits, subject, query, type)
{
  # test (up to 'max') whether the 5' stretches of the subject
  # and the query are identical
  o <- apply(hits, 1, function(x){
    if( (x[8] == "+" & type == '5p') ){
      s <- try( stringi::stri_sub(subject[x[1]], 1, 1:(as.numeric(x[6]) - 1)), silent = TRUE )
    } else if( (x[8] == "-" & type == '3p') ){
      s <- try( reverseComplement(stringi::stri_sub(subject[x[1]], 1, 1:(as.numeric(x[6]) - 1))), silent = TRUE )
    } else if((x[8] == "-" & type == '5p') ){
      s <- try( reverseComplement(stringi::stri_sub(subject[x[1]], as.numeric(x[7]), (as.numeric(x[7])):nchar(subject[x[1]]))),
                silent = TRUE )
    } else if( (x[8] == "+" & type == '3p') ){
      s <- try( stringi::stri_sub(subject[x[1]], as.numeric(x[7]), (as.numeric(x[7])):nchar(subject[x[1]])), silent = TRUE )
    }
    if(class(s) == "try-error") return(NA)
    query <- gsub(".", "", query, fixed = TRUE)
    if (type == '5p') {
      q <- try( stringi::stri_sub( query, 1, 1:(as.numeric(x[4]) - 1) ), silent = TRUE )
    } else if(type == '3p') {
      q <- try( stringi::stri_sub( query, as.numeric(x[5]) + 1, (as.numeric(x[5]) + 1):nchar(query) ), silent = TRUE )
    }
    if( class(q) == "try-error" ) return(NA)
    suppressWarnings( sum(stringi::stri_cmp_eq(tolower(q), tolower(s))) )
  })
  names(o) <- hits$subject
  o
}

