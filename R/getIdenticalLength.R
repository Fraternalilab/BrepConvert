#' Function to calculate the length of identical sequence stretch 5' or 3' of given boundaries of sequence stretches from two longer sequences.
#'
#' @param hits data.frame storing potential hits (donor pseudogenes) that explain a given gene conversion event. Each row corresponds to one definition using a given pseudogene. The boundaries of the pseudogene sequence and the observed repertorie sequence are given.
#' @param query character, the full-length V gene sequence of the observed repertoire sequence for which gene conversion events are annotated.
#' @param subject \code{DNAStringSet} object containing a set of pseudogene DNA sequences.
#' @param type character, either '5p' or '3p' denoting the direction for which the identical sequence stretches are considered.
#' 
#' @importFrom stringi stri_sub stri_cmp_eq
#' @importFrom Biostrings reverseComplement nchar
#'
getIdenticalLength <- function(hits, subject, query, type)
{
  # test (up to 'max') whether the 5' stretches of the subject
  # and the query are identical
  o <- apply(hits, 1, function(x){
    if( (x[8] == "+" & type == '5p') ){
      s <- try( rev(stringi::stri_sub(subject[x[1]], 1:(as.numeric(x[6]) - 1), (as.numeric(x[6]) - 1))), silent = TRUE )
    } else if( (x[8] == "-" & type == '3p') ){
      s <- try( rev(Biostrings::reverseComplement(stringi::stri_sub(subject[x[1]], 1:(as.numeric(x[6]) - 1), (as.numeric(x[6]) - 1)))), silent = TRUE )
    } else if((x[8] == "-" & type == '5p') ){
      s <- try( Biostrings::reverseComplement(stringi::stri_sub(subject[x[1]], as.numeric(x[7]) + 1, (as.numeric(x[7]) + 1):Biostrings::nchar(subject[x[1]]))),
                silent = TRUE )
    } else if( (x[8] == "+" & type == '3p') ){
      s <- try( stringi::stri_sub(subject[x[1]], as.numeric(x[7]) + 1, (as.numeric(x[7]) + 1):Biostrings::nchar(subject[x[1]])), silent = TRUE )
    }
    if(class(s) == "try-error") return(NA)
    query <- gsub(".", "", query, fixed = TRUE)
    if (type == '5p') {
      q <- try( rev(stringi::stri_sub( query, 1:(as.numeric(x[4]) - 1), (as.numeric(x[4]) - 1) )), silent = TRUE )
    } else if(type == '3p') {
      q <- try( stringi::stri_sub( query, as.numeric(x[5]) + 1, (as.numeric(x[5]) + 1):nchar(query) ), silent = TRUE )
    }
    if( class(q) == "try-error" ) return(NA)
    suppressWarnings( sum(stringi::stri_cmp_eq(tolower(q), tolower(s))) )
  })
  names(o) <- hits$subject
  o
}

