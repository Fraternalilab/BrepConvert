## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(BrepConvert)

dataset <- system.file( "extdata/pmid3100050_vquest.tsv", 
                        package = "BrepConvert" )
dataset <- read.table( dataset, sep = "\t", stringsAsFactors = FALSE, header = TRUE )
dataset[, c("sequence_id", "locus", "v_call", "j_call")]

# 'sequence_alignment' contain the IMGT-gapped DNA sequence of the V gene 
# for each observation
print( head( dataset$v_sequence_alignment, 3 ) )

# assigning it to a new variable and name this vector with the sequence_id
repertoire <- dataset$v_sequence_alignment
names( repertoire ) <- dataset$sequence_id
print( head( repertoire , 3) )

## -----------------------------------------------------------------------------
# FASTA files containing pseudogene and functional alleles of the Chicken IGLV
# locus are shipped with the package.
functional_IGLV <- system.file("extdata/IMGT_Chicken_IGLV_F.fasta",
                               package = "BrepConvert")
pseudogene_IGLV <- system.file("extdata/IMGT_Chicken_IGLV_P.fasta",
                               package = "BrepConvert")

# An executable of the BLAT program is also shipped with the package.
blat <- system.file("exe/blat", package = "BrepConvert")
# NOTE: This works for Linux OS only. For other OS, please download/compile
# the executable and indicate the filepath like so:
# blat <- "/home/abc/Documents/blat/blat"

annotation <- batchConvertAnalysis(
  functional = functional_IGLV,
  pseudogene = pseudogene_IGLV,
  repertoire = repertoire, # notice here the vector is passed, NOT the entire table!
  blat_exec = blat
)

# look at the annotation for sequence '3W-3'
print( annotation[ annotation$SeqID == "3W-3", ] )


## ---- fig.width=7, fig.height=6-----------------------------------------------
plotEvents( annotation = annotation, repertoire = repertoire,
            show.sequence.names = TRUE, # recommended FALSE if you are visualising a lot of sequences!
            thickness = c(4, 1), # this controls the thickness of lines drawn
            highlight.genes = c("IGLV1-2", "IGLV1-16", "IGLV1-22" )) # this is optionally; these donors will be highlighted with different colours

## ---- fig.width=7, fig.height=7-----------------------------------------------
# first add an extra column 'sequence_type' indicating whether
# a sequence is from the 3W or 18D group
annotation$sequence_type <- ifelse(grepl("^3W-", annotation$SeqID), "3W", "18D")
plotCoverage(annotation, repertoire, 'sequence_type', 
             type = 'coverage', # alternatively 'start' considers only the start sites
             definition = 'both') # or 'narrow' / 'broad' if you only want either definition of boundaries

