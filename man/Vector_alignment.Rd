% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/array_structrue.R
\name{Vector_alignment}
\alias{Vector_alignment}
\alias{summary_SV}
\title{Dynamic programming, a little slow but can used for multi-characteristics.}
\usage{
Vector_alignment(seq1, seq2, mat = NULL, gap = 0, print_score_mat = FALSE)

summary_SV(align_res)
}
\arguments{
\item{seq1}{vector}

\item{seq2}{vector}

\item{mat}{substitution matrix}

\item{gap}{gap penalty}

\item{print_score_mat}{print score matrix}

\item{align_res}{result from Vector_alignment}
}
\value{
alignment result

list
}
\description{
Dynamic programming, a little slow but can used for multi-characteristics.

Summary SV of an alignment
}
\examples{
align_res <- Vector_alignment(c("1" = "as", "2" = "bb", "3" = "cc"), c("bb", "cc", "as", "as"))
align_res
summary_SV(align_res)
Vector_alignment("ATCGTACG", "ATCGTAGC", print_score_mat = TRUE)
}
