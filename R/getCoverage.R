#' Count gene conversion events which cover/start at every DNA position along the germline
#'
#' @param db data.frame, output from \code{batchConvertAnalysis}.
#' @param column1 optional, character corresponding to a column name in \code{db}. If supplied, calculation of coverage will be considered separately for every unique value in \code{column1} of the data frame \code{tb} so that histograms can be plotted separately.
#' @param column2 optional, character corresponding to a column name in \code{db}. If supplied, calculation of coverage will be considered separately for every unique value in \code{column2} of the data frame \code{tb} so that histograms can be plotted separately.
#' @param start_col character, column in \code{db} which gives the start position.
#' @param end_col character, column in \code{db} which gives the end position. Not used when \code{type == 'start'}.
#' @param type character, either 'coverage' (i.e. how often is a given DNA position included within gene conversion events?) or 'start' (i.e. how often do gene conversion events start at a given DNA position?).
#'
#' @description This function counts how often gene conversion events (partitioned by different combinations of values in \code{column1} and \code{column2} of \code{db}) cover/start at every DNA position along the V gene.
#'
#' @return A data.frame with three columns:
#' \describe{
#'   \item{position}{numeric, DNA position}
#'   \item{value}{numeric, if \code{type =='coverage'} the proportion of gene conversion events which cover the given \code{position}. If \code{type == 'start'}, the proportion of events which begin at the given \code{position}.}
#'   \item{region}{character, either 'FWR' (framework) or 'CDR'. Follows the standard IMGT numbering; the region is inferred from the \code{position} column.}
#' }
#'
getCoverage <- function(db, column1, column2, start_col,
                        end_col = NULL, type = c('coverage', 'start'))
{
  if ( ! is.data.frame(db) )
    stop( "'db' must be a data.frame, output from batchConvertAnalysis." )
  if( ! type %in% c('coverage', 'start'))
    stop("'type' must be either 'coverage' or 'start'.")
  if( is.null( column1 ) ){
    c1all <- 1
  } else {
    if( ! column1 %in% colnames( db ))
      stop("'column1' must be a column named in 'db'.")
    c1all <- unique(as.character(db[, column1]))
    c1all <- c1all[ !is.null(c1all) ]
    c1all <- c1all[ !is.na(c1all) ]
  }
  if( is.null( column2 ) ){
    c2all <- 1
  } else {
    if( ! column2 %in% colnames( db ))
      stop("'column1' must be a column named in 'db'.")
    c2all <- unique(as.character(db[, column2]))
    c2all <- c2all[ !is.null(c2all) ]
  c2all <- c2all[ !is.na(c2all) ]
  }
  if( !is.character( start_col ))
    stop("'start_col' should be a character string.")
  if( !is.character( end_col ))
    stop("'end_col' should be a character string.")
  if( ! start_col %in% colnames( db ) )
    stop("'start_col' should be a column named in 'db'. If you are indicating the 'soft_start' i.e. start site padded with flanking sequence, please make sure you are supplying the output of adjustBroadPos to this function.")
  if( ! end_col %in% colnames( db ) )
    stop("'end_col' should be a column named in 'db'. If you are indicating the 'soft_end' i.e. start site padded with flanking sequence, please make sure you are supplying the output of adjustBroadPos to this function.")
  covg <- do.call("rbind", lapply(c1all, function(c1){
    do.call("rbind", lapply(c2all, function(c2){
      if( length(c1all) == 1 && length(c2all) == 1 && c1all == 1 && c2all == 1){
        tb2 <- db
      } else if( length(c1all) == 1 && c1all == 1 ){
        tb2 <- db[which(db[, column2] == c2), ]
      } else if( length(c2all) == 1 && c2all == 1 ){
        tb2 <- db[which(db[, column1] == c1), ]
      }
      ranges <- IRanges::IRanges(start = tb2[, start_col], end = tb2[, end_col])
      if(type == "coverage"){
        covg <- as.data.frame(IRanges::coverage(ranges))
        covg[, 1] <- covg[, 1] / nrow(tb2)
        covg$position <- 1:nrow(covg)
      } else if(type == "start"){
        covg <- plyr::ddply(tb2, start_col, nrow)
        covg[, 2] <- covg[, 2] / nrow(tb2)
        colnames(covg) <- c("position", "value")
      }
      covg[, column2] <- c2
      covg[, column1] <- c1
      covg
    }))
  }))
  # cap at position 330 (ie AA position 110 - roughly the end of V genes in BCRs)
  covg <- covg[ which( covg$position <= 330 ), ]
  covg$region <- sapply(covg$position, function(x){
    # IMGT defs
    if( x %in% c(1:78, 115:165, 196:312) ) return("FW")
    if( x %in% c(79:114, 166:195, 313:330) ) return("CDR")
  })
  covg
}
