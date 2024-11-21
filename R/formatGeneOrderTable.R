#' Format Gene Order Table from IMGT into a data frame
#'
#' @param html_table raw html code output from \code{getGeneOrderTable}.
#'
#' @description
#' This function parse the HTML code scraped using the \code{geneGeneOrderTable} function into a data frame.
#' Due to the fact that IMGT formats these tables in different species slightly differently (they have different number of genome assemblies, haplotypes etc.), the exact number of columns (and their names) depend on the queried species and locus.
#'
#' @return A data frame with separate columns indicating immunoglobulin gene names and gene order (as an integer, with 1 being the most 5' gene in the locus).
#'
#' @importFrom rvest html_text html_nodes

# Function to format the HTML table into a data frame
formatGeneOrderTable <- function(html_table) {
  # Extract all rows from the table
  rows <- rvest::html_nodes(html_table, "tr")

  # Initialize an empty list to store the rows
  table_list <- list()

  # Track the maximum number of columns found in any row
  max_columns <- 0

  # Loop through each row and process its cells
  for (i in seq_along(rows)) {
    # Extract the individual cells in the row (both <th> and <td>)
    cells <- rvest::html_nodes(rows[[i]], xpath = ".//th | .//td")

    # Extract the text content from each cell
    row_values <- rvest::html_text(cells)
    # special case is for the human IGH table on IMGT
    # strange formatting. Here set an arbitrary cutoff that if
    # row has only 1 column (in this case just a \t)
    # or with many (here catch this with >1000) I am ignoring them
    if(length(row_values) == 1){
      next
    }
    if(length(row_values) > 1000){
      next
    }

    # Track the maximum number of columns seen
    max_columns <- max(max_columns, length(row_values))

    # Store the row in the table list
    table_list[[i]] <- row_values
  }

  # Now, ensure that each row has the same number of columns (fill with NA if needed)
  table_list <- lapply(table_list, function(row) {
    # If a row has fewer columns than the maximum, fill it with NA at the end
    if (length(row) < max_columns) {
      row <- c(row, rep(NA, max_columns - length(row)))
    }
    return(row)
  })

  # Convert the list of rows into a data frame
  table_df <- do.call(rbind, lapply(table_list, function(x) as.data.frame(t(x), stringsAsFactors = FALSE)))

  # Clean the column names (using the first row as headers)
  colnames(table_df) <- make.names(unlist(table_list[[1]]))  # Use the first row as headers
  table_df <- table_df[-1, ]  # Remove the header row

  table_df <- table_df[apply(table_df, 1, function(x) sum(is.na(x))) < ncol(table_df), ]

  # remove space from strings
  for(i in 1:ncol(table_df)){
    if(class(table_df[, i]) == "character"){
      table_df[, i] <- gsub(" ", "", table_df[, i])
    }
  }

  # Return the cleaned data frame
  return(table_df)
}


