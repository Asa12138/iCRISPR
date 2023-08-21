# Load the value of the option on package startup
.onAttach <- function(libname, pkgname) {
  if(!dir.exists(tools::R_user_dir("iCRISPR")))dir.create(tools::R_user_dir("iCRISPR"),recursive = TRUE)
  refresh_config()
}

# read user config first.
refresh_config <- function() {
  default_options <- readRDS(file = system.file("config", package = "iCRISPR"))
  file_path=file.path(tools::R_user_dir("iCRISPR"),"config")
  if (file.exists(file_path)) {
    options_to_load <- readRDS(file = file_path)
    if(length(options_to_load)==0)options_to_load=NULL
    options_to_load=pcutils::update_param(default_options,options_to_load)
  }
  else {
    options_to_load=default_options
  }
  # set options
  options("iCRISPR_config"=options_to_load)
}

#' Show config
#'
#' @return config
#' @export
show_config=function(){
  config=getOption("iCRISPR_config")
  return(config)
}

#' Set config
#'
#' @param item item
#' @param value value
#'
#' @return No value
#' @export
#'
set_config <- function(item,value) {
  refresh_config()
  config=getOption("iCRISPR_config")
  if(is.null(value))config=config[-which(names(config)==item)]
  else config=pcutils::update_param(config,setNames(list(value),item))
  saveRDS(config, file = file.path(tools::R_user_dir("iCRISPR"),"config"))
  options("iCRISPR_config"=config)
  message("Set sucessfully!")
}


#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
#' @param lhs A value or the magrittr placeholder.
#' @param rhs A function call using the magrittr semantics.
#' @return The result of calling `rhs(lhs)`.
