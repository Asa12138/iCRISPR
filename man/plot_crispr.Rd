% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CRISPRCasfinder_downstream.R
\name{plot_crispr}
\alias{plot_crispr}
\title{Plot a crispr-cas system}
\usage{
plot_crispr(crispr, genome = NULL, contig = NULL)
}
\arguments{
\item{crispr}{crispr result from \code{pre_CCF_res()}}

\item{genome}{the genome name}

\item{contig}{the contig name}
}
\value{
ggplot
}
\description{
Plot a crispr-cas system
}
\examples{
data(crispr)
plot_crispr(crispr,genome="MAG_test",contig="AAB-S01R1_k55_9399631_flag=0_multi=9.8751_len=26518")
}