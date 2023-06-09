% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Seqinfo.R
\name{summary_seq}
\alias{summary_seq}
\title{Summarize sequence information}
\usage{
summary_seq(sequence, max_k = 7)
}
\arguments{
\item{sequence}{The input DNA sequence.}

\item{max_k}{The maximum k value for kmers.}
}
\description{
Summarize sequence information
}
\examples{
summary_seq("CCTGAACCTATGCCGTCCACCTTGCGTTGCCTT")
apply(random_seq(10),1,\(i)summary_seq(i[2]))

}