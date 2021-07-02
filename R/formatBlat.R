#' Function to format BLAT output
#'
#' @param tb data.frame, output from the \code{blat} function.
#' @param adjustment integer, the number of nucleotides which flank the gene conversion events (included for BLAT just to facilitate the alignment to look for seeds). This needs to be taken into account to re-adjust the actual start and end point of conversion events. (default: 6)
#' 
#' @examples
#' \dontrun{ }
#'
#' @importFrom stringr str_extract
#' 
formatBlat <- function(tb, adjustment = 6)
{
  tb[, "gene"] <- gsub("|", "", stringr::str_extract(tb[, "gene"], "IG.V.*\\||IG.V.*$"),
                       fixed = TRUE)
  tb[, "g_start"] <- sapply(tb$qname, function(x){
    y <- unlist(strsplit(x, split = "/"))
    as.numeric( unlist(strsplit(y[length(y)], split = "-"))[1] )
  })
  tb[, "q_start"] <- tb[, "q_start"] + tb[, "g_start"] - 1 - adjustment
  tb[, "q_end"] <- tb[, "q_end"] + tb[, "g_start"] - 1 - adjustment
  tb[, "qname"] <- sapply(tb$qname, function(x){
    y <- unlist(strsplit(x, split = "/"))
    paste( y[-length(y)], collapse = "/" )
  })
  tb <- tb[, c("qname", "gene", "mismatch", "gaps", "q_start", "q_end",
               "s_start", "s_end")]
  if( is.data.frame(tb) ){
    if( nrow(tb) > 0 ){
      colnames(tb) <- c("query", "subject", "mismatches", "gaps_open", "q_start", "q_end",
                        "s_start", "s_end")
      return(tb)
    } else return(data.frame())
  } else return(data.frame())
}
