% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generateSeq.R
\name{generateSeq}
\alias{generateSeq}
\title{generateSeq
generate the isoform sequence from the isoform list}
\usage{
generateSeq(genome, isoform)
}
\arguments{
\item{genome}{DNAStringSet indicating the genome reference}

\item{isoform}{a list of data frame indicating the exon locations in isoforms}
}
\value{
a list comprsing the following items 
\itemize{
 \item {rna} character value indicating the RNA sequence of isoforms
 \item {protein} character value indicating the amino acid sequence of 
 isoforms
 \item {translated}  boolean value indicating whether the isoform can be 
 translated 
}
}
\description{
generateSeq
generate the isoform sequence from the isoform list
}
