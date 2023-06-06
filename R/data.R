#'@title One crispr object result from MAG_test.
#'@description a list contains four dataframe.
#'
#'@docType data
#'@name crispr
#'@usage crispr
#'@format contains CRISPR, Cas, Array and Spacer
#'\describe{
#' \item{CRISPR}{CRISPR found in the whole genome (with multi-contigs), CRISPR_id is the key}
#' \item{Cas}{Cas system found in the whole genome, Cas_id is the key, one Cas_id contains multi cas protein}
#' \item{Array}{the array structure of each CRISPR, "LeftFLANK","CRISPRdr","CRISPRspacer","RightFLANK"}
#' \item{Spacer}{all spacer come from each CRISPR, Spacer_id is the key}
#'}
#'
NULL
