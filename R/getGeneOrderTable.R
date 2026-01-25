#' Fetch order of immunoglobulin genes at the genomic locus
#'
#' @param species string, common name of the species of interest (see 'Description' for further details).
#' @param locus string, immunoglobulin locus (IGH, IGK, IGL) (see 'Description' for further details).
#'
#' @description
#' This function fetch from the IMGT website on the order of immunoglobulin genes at a given genomic locus for a user-specified species.
#' The list of species can be viewed here: https://www.imgt.org/IMGTrepertoire/LocusGenes/ (under "Locus gene order"). The common names of species were expected in this function.
#'
#' @return A data frame with separate columns indicating immunoglobulin gene names and gene order (as an integer, with 1 being the most 5' gene in the locus).
#'
#' @importFrom stringr str_to_lower
#' @importFrom rvest read_html html_nodes
#' @export getGeneOrderTable
getGeneOrderTable <- function(species, locus) {
  if(! locus %in% c("IGH", "IGK", "IGL")){
    message("locus must be any one of 'IGH', 'IGL' or 'IGK'.")
    return(NULL)
  }
  url <- paste0("https://www.imgt.org/IMGTrepertoire/index.php?section=LocusGenes&repertoire=GeneOrder&species=",
                gsub(" ", "_", stringr::str_to_lower(species)), "&group=", locus)

  # Read the HTML content from the URL
  webpage <- try(rvest::read_html(url), silent = TRUE)
  if(any(class(webpage) == "try-error")){
    message("Cannot read data. Are you sure the supplied species is correct?")
    return(NULL)
  }
  # Find all tables on the webpage
  tables <- rvest::html_nodes(webpage, "table")

  # If there are no tables, return NULL
  if (length(tables) == 0) {
    message("No tables found on the webpage.")
    return(NULL)
  }

  # Extract all tables and convert it to a data frame
  out <- lapply(tables, formatGeneOrderTable)

  if(length(out) == 1) return(out[[1]])
  if(species == "human" & locus == "IGH"){
    # for some reason there are two tables but they have the same information
    # only keep the first
    out <- out[[1]]
    out <- out[, 1:4] # only the first 4 columns are useful
    # first row has legit column names
    colnames(out) <- out[1, ]
    return(out[-1:-4, ]) # first 3 rows don't have relevant info
  }
  return(out)
}


