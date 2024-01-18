#' @title One crispr object result from MAG_test.
#' @description a list contains four dataframe.
#'
#' @docType data
#' @name crispr
#' @format contains CRISPR, Cas, Array and Spacer
#' \describe{
#' \item{CRISPR}{CRISPR found in the whole genome (with multi-contigs), CRISPR_id is the key}
#' \item{Cas}{Cas system found in the whole genome, Cas_id is the key, one Cas_id contains multi cas protein}
#' \item{Array}{the array structure of each CRISPR, "LeftFLANK","CRISPRdr","CRISPRspacer","RightFLANK"}
#' \item{Spacer}{all spacer come from each CRISPR, Spacer_id is the key}
#' }
#'
NULL

#' @title One multi_crispr object contains lots of crispr object.
#' @description a list contains lots of crispr object (50 genomes).
#'
#' @docType data
#' @name multi_crispr
#'
NULL

#' @title One multi_crispr object contains lots of crispr object.
#' @description a list contains lots of crispr object (50 genomes).
#'
#' @docType data
#' @name multi_crispr2
#'
NULL

#' @title CRISPR-Cas system type information.
#'
#' @docType data
#' @name cas_type_df
#'
NULL

#' @title CRISPR-Cas system type tree.
#'
#' @docType data
#' @name cas_tree
#'
NULL

#' @title Prokaryotic CRISPR source-target network
#'
#' @docType data
#' @name pro_net
#'
NULL
