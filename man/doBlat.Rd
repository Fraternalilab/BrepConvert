% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/doBlat.R
\name{doBlat}
\alias{doBlat}
\title{Extract sequence stretches to be analysed using BLAT.}
\usage{
doBlat(seqname, repertoire, functional, gapwidth = 3, flank = 6)
}
\arguments{
\item{seqname}{character, name of sequence in \code{repertoire}.}

\item{repertoire}{a named vector of full-length sequences to be analysed. The \code{names} attribute of this vector stores the sequence identifiers.}

\item{functional}{\code{DNAStringSet} object storing DNA sequences of the functional V gene alleles.}

\item{gapwidth}{integer, the distance between separate mismatches to be grouped together as a 'gap' to be analysed as candidate gene conversion events. (default: 3)}

\item{flank}{integer, the number of nucleotides to be included in the 'gap' sequences to be searched using BLAT to annotate for candidate gene conversion events. (default: 6)}
}
\description{
Extract sequence stretches to be analysed using BLAT.
}
