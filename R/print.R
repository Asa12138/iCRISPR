#' Print crispr object
#'
#' @method print crispr
#' @param ... add
#' @param x crispr object
#'
#' @exportS3Method
#' @return No value
#' @examples
#' data(crispr)
#' print(crispr)
print.crispr <- function(x, ...) {
  info <- attributes(x)$basic_info
  pcutils::dabiao("Genome name: ", info$genome_name)
  cat("With", info$n_crispr, "CRISPR systems(", sum(x$CRISPR$evidence_level == 4), "evidence level=4)\n")
  cat("With", info$n_cas, "CAS systems\n")
  cat("With", info$n_spacer, "spacers\n")
}


#' Print multi_crispr object
#'
#' @method print multi_crispr
#' @param ... add
#' @param x multi_crispr object
#'
#' @exportS3Method
#' @return No value
print.multi_crispr <- function(x, ...) {
  pcutils::dabiao("a list with total ", length(x), " genomes", print = TRUE)
}


#' Print Array object
#'
#' @method print Array
#' @param ... add
#' @param x Array object
#'
#' @exportS3Method
#' @return No value
print.Array <- function(x, ...) {
  pcutils::dabiao("CRISPR_id: ", x$CRISPR_id, print = TRUE)
  cat("With", x$spacer_number, "spacers; Evidence level=", x$evidence_level, "\n")
  cat("With consensus_repeat:", x$consensus_repeat, "\n")
}
