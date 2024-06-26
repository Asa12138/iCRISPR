% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CRISPRCasfinder_downstream.R
\name{plot_crispr}
\alias{plot_crispr}
\alias{xlim_crispr}
\title{Plot a crispr-cas system}
\usage{
plot_crispr(crispr, contig = NULL, genome = NULL, array = TRUE, cas = TRUE)

xlim_crispr(...)
}
\arguments{
\item{crispr}{crispr result from \code{pre_CCF_res()}}

\item{contig}{the contig name}

\item{genome}{the genome name}

\item{array}{plot the array?}

\item{cas}{plot the cas?}

\item{...}{numeric}
}
\value{
ggplot

ggproto object
}
\description{
Plot a crispr-cas system

xlim for plot_crispr
}
\examples{
data(crispr)
(p <- plot_crispr(crispr,
  genome = "MAG_test",
  contig = "AAB-S01R1_k55_9399631_flag=0_multi=9.8751_len=26518"
))
p + xlim_crispr(24000, 27000)
}
