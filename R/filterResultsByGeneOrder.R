#' Visualise locations of gene conversion events
#'
#' @param results data.frame, output from \code{batchConvertAnalysis}.
#' @param gene_order A named vector of numeric/integers with the order of germline genes in the genome. Can be obtained from IMGT using the \code{getGeneOrderTable} function in BrepConvert, and further generating this vector from the gene order table (see package vignette for example). Each element should be named by the corresponding germline gene name.
#' @param functional_alelle_col Column name in \code{results} denoting the annotated functional allele for the sequence. (Default: 'allele')
#' @param donor_allele_col Column name in \code{results} with semicolon-delimited list of annotated pseudogenes for the detected gene conversion event.
#'
#' @description This function takes output from the \code{batchConvertAnalysis} function, and add a new column "filtered_pseudogenes" which is a pruned version of the \code{donor_allele_col} removing pseudogenes which are physically implausible (because they are 3' of the acceptor gene - they would have been removed from the genome during V(D)J recombination).
#'
#' @return The dataframe in \code{results} with a new column named 'filtered_pseudogenes' appended. This column contains the pruned list of possible genes as donors for the observed gene conversion event, after accounting for the physical order of immunoglobulin gene on the chromosome.
#'
#' @export filterResultsByGeneOrder
filterResultsByGeneOrder <- function(results, gene_order,
                                     functional_alelle_col = "allele", donor_allele_col = "gene")
{
  if( is.null(names(gene_order) ))
    stop("'gene_order' should be a vector of integers with names corresponding to each germline gene.")
  results <- split(results, f = results[, functional_alelle_col])
  results <- lapply(results, function(tb){
    functional_allele <- unlist(strsplit(tb[1, functional_alelle_col], split = "*", fixed = TRUE))[1]
    if(! unlist(strsplit(functional_allele, split = "*", fixed = TRUE))[1] %in% names(gene_order)){
      warning(paste0(functional_allele, " is not included in 'gene_order'. Not filtering any results for this gene."))
      return(NULL)
    } else {
      j = gene_order[functional_allele]
      donor_genes_to_keep <- names(gene_order[gene_order <= j])
      tb$filtered_pseudogenes <- sapply(tb[, donor_allele_column], function(x){
        if(is.na(x)) return(x)
        all_genes <- unlist(strsplit(x, split = ";"))
        all_genes <- unique(sapply(all_genes, function(y){
          unlist(strsplit(y, split = "*", fixed = TRUE))[1]
        }))
        return(paste(all_genes[grepl(paste(donor_genes_to_keep, collapse = "|"), all_genes)],
                     collapse = ";"))
      })
      return(tb)
    }
  })
  do.call("rbind", results)
}
