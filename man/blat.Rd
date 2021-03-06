% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/blat.R
\name{blat}
\alias{blat}
\title{Function to run BLAT.}
\usage{
blat(sequenceTb, database, blat_exec, min_score = 20)
}
\arguments{
\item{sequenceTb}{data.frame, each row corresponding to one sequence to be queried using BLAT. The data frame consists of the following two columns:
\describe{
  \item{qname}{character, name/identifier of the query sequence.}
  \item{seq}{character, (nucleotide) sequence to be queried using BLAT.}
}}

\item{database}{character, filepath to the FASTA file storing the sequence database to be searched against using BLAT.}

\item{blat_exec}{character, filepath to the executable of the BLAT program.}

\item{min_score}{integer, minimum score to trigger a BLAT alignment, i.e. the \code{minScore} parameter in the command-line BLAT program (default: 20).}
}
\description{
This function invokes BLAT from the command line to search sequences given in \code{sequenceTb}, against FASTA sequences stored in the file \code{database}. The parameters used are: \code{-out=blast8 -tileSize=6 -stepSize=1}, and \code{minScore} given by the parameter \code{min_score} of this function. The results are written to a temporary file and then read in as a data frame.
}
