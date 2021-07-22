#' Function to run BLAT.
#'
#' @param sequenceTb data.frame, each row corresponding to one sequence to be queried using BLAT. The data frame consists of the following two columns:
#'   \describe{
#'     \item{qname}{character, name/identifier of the query sequence.}
#'     \item{seq}{character, (nucleotide) sequence to be queried using BLAT.}
#'   }
#' @param database character, filepath to the FASTA file storing the sequence database to be searched against using BLAT.
#' @param blat_exec character, filepath to the executable of the BLAT program.
#' @param min_score integer, minimum score to trigger a BLAT alignment, i.e. the \code{minScore} parameter in the command-line BLAT program (default: 20).
#'
#' @description This function invokes BLAT from the command line to search sequences given in \code{sequenceTb}, against FASTA sequences stored in the file \code{database}. The parameters used are: \code{-out=blast8 -tileSize=6 -stepSize=1}, and \code{minScore} given by the parameter \code{min_score} of this function. The results are written to a temporary file and then read in as a data frame.
#'
#' @importFrom utils read.table
#'
blat <- function(sequenceTb, database, blat_exec, min_score = 20)
{
  if(!file.exists(blat_exec)){
    stop("BLAT executable does not exist. Are you sure you are passing the right filepath?")
  }
  checkBlat <- suppressWarnings( system(blat_exec, intern = TRUE) )
  if( !grepl("blat - Standalone BLAT", checkBlat[1], fixed = TRUE) ){
    stop("blat_exec does not point to the BLAT executable. Are you sure you are passing the right filepath?")
  }
  if(!is.integer(as.integer(min_score))){
    stop("min_score should be an integer.")
  }
  if(!is.data.frame(sequenceTb)){
    stop("sequenceTb should be a data frame.")
  }
  if(!all(c("qname", "seq") %in% colnames(sequenceTb))){
    stop("Check that sequenceTb contains all of the following columns: qname, seq")
  }

  # check system
  if (.Platform$OS.type == "windows") {
    invoke <- shell
  } else {
    invoke <- system
  }

  sequence_fasta <- tempfile()
  for(i in 1:nrow(sequenceTb)){
    write( c(paste0(">", sequenceTb[i, "qname"]), as.character(sequenceTb[i, "seq"])) ,
           file = sequence_fasta, append = TRUE )
  }
  #minimum 6-nt match to trigger an alignment - minimum for nt sequences
  #minScore=20 to catch exons (identical to web-browser version of blat)
  random_n <- sample(1:10000, 1)
  min_score <- as.integer(min_score)
  cmd <- paste0(blat_exec, ' ', database, " ", sequence_fasta,
                " blat_output", random_n, ".blast -out=blast8 -tileSize=6 -stepSize=1 -minScore=", min_score)
  do.call(invoke, list("command" = cmd, "intern" = FALSE , "ignore.stdout" = TRUE, "ignore.stderr" = TRUE))
  result_table <- try( read.table(paste0("blat_output", random_n, ".blast"), stringsAsFactors = F),
                       silent = TRUE)
  file.remove(sequence_fasta)
  if( class(result_table) != "try-error" ){
    out <- result_table[, c(1, 7, 8, 9, 10, 5, 6, 2)] #, ncol(result_table) - 1, ncol(result_table))]
    colnames(out) <- c("qname", "q_start", "q_end", "s_start", "s_end",
                       "mismatch", "gaps", "gene")#, "CB", "UB")
    out[, "dir"] <- apply(result_table, MARGIN = 1, function(x){
      if( ((as.numeric(x[7]) < as.numeric(x[8])) & (as.numeric(x[9]) > as.numeric(x[10]))) |
          ((as.numeric(x[7]) > as.numeric(x[8])) & (as.numeric(x[9]) < as.numeric(x[10]))) ){
        return("-")
      } else return("+")
    })
    file.remove(paste0("blat_output", random_n, ".blast"))
    return( out[, c("qname", "q_start", "q_end", "s_start", "s_end", "mismatch",
                    "gaps", "dir", "gene")] )
  } else return(NULL)
}
