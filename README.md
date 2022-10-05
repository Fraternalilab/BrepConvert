# BrepConvert
R package to analyse gene conversion events in immunoglobulin sequences

## Install

in R:

```
require(devtools)
install_github('Fraternalilab/BrepConvert', ref = 'main')
```

## How to use

This [vignette](http://htmlpreview.github.io/?https://github.com/Fraternalilab/BrepConvert/blob/main/vignettes/BrepConvert.html) contains examples on how to use `BrepConvert.`

The analysis performed on Mallby et al. was performed as follows: IMGT-gapped V gene sequences were splitted into batches of 3000 sequences, each read into R as a named vector. One of these files is uploaded here in the directory `example`.

The following shows the run-code to generate the result of this file via BrepConvert:

```
library(Biostrings)
library(BrepConvert)

sequences <- readDNAStringSet("MallabyEtAl_example_set.fasta")
functional_IGHV <- "IGHV_functional_gapped.fasta" # Downloads 21/08/2020
pseudogene_IGHV <- "IGHV_Pseudogenes_gapped.fasta"

blat <- system.file("exe/blat", package = "BrepConvert")

annotation <- batchConvertAnalysis(
  functional = functional_IGHV,
  pseudogene = pseudogene_IGHV,
  repertoire = sequences,
  blat_exec = blat
)
```

The corresponding output file (in CSV format) is also uploaded into the `example` directory in this repository.

# Cite

Mallby J, Ng JCF, Stewart AS, Sinclair E, Dunn-Walters DK, Hershberg U. Diversification of Antibodies by Gene Conversion in the Domestic Chicken (*Gallus gallus domesticus*). *Under Review*.
