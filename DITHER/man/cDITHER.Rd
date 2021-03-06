% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cDITHER.R
\name{cDITHER}
\alias{cDITHER}
\title{cDITHER}
\usage{
cDITHER(input_data_cna)
}
\arguments{
\item{input_data_cna}{A dataframe or matrix which at least includes two columns ("Sample" and "Segment_Mean").}
}
\value{
A dataframe with 2 columns:
\item{Sample}{Tumor sample ID.}
\item{cDITHER_score}{The CNA entropy-based ITH score in tumors.}
}
\description{
\code{cDITHER} Calculate ITH based on the entropy of CNA profiles in tumors.
}
\details{
This function calculates the ITH score in a bulk tumor based on the entropy of its CNA profiles.
}
\examples{
path = system.file("extdata", "example_cDITHER.txt", package = "DITHER", mustWork = TRUE)
input_data_cna = read.table(path, stringsAsFactors = FALSE, sep = "\\t", header = TRUE, quote = "")
cDITHER(input_data_cna)
}
\author{
Lin Li <cpu_lilin@126.com>, Xiaosheng Wang <xiaosheng.wang@cpu.edu.cn>
}
