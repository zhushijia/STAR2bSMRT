% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/starShort.R
\name{starShort}
\alias{starShort}
\title{starShort
run the STARshort mapping}
\usage{
starShort(genomeDir, SR1, SR2, outputDir, SJ = NULL)
}
\arguments{
\item{genomeDir}{character value indicating the location of STAR genome indexes}

\item{SR1}{character value indicating the short read file in fastq format: 
single-end or paired-end R1}

\item{SR2}{character value indicating the short read file in fastq format: 
paired-end R2}

\item{outputDir}{character value indicating the directory of output}

\item{SJ}{character value for 2pass STAR mapping}
}
\description{
starShort
run the STARshort mapping
}
