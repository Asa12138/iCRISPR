# read user config first.
refresh_config=function(){
  default_config=read.csv(system.file("config", package = "iCRISPR"),header =TRUE,quote = "")
  default_config=setNames(default_config$value,default_config$item)
  if(file.exists(file.path(tools::R_user_dir("iCRISPR"),"config"))){
    config=read.csv(file.path(tools::R_user_dir("iCRISPR"),"config"),header =TRUE,quote = "")
    if(nrow(config)>0){
      config=setNames(config$value,config$item)
    }
    else config=NULL
    config=pcutils::update_param(default_config,config)
  } else {
    config=default_config
  }
  config
}

#' Show config
#'
#' @return config
#' @export
show_config=function(){
  config=refresh_config()
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
set_config=function(item,value){
  config=refresh_config()
  if(is.null(value))config=config[-which(names(config)==item)]
  else {
    config[item]=value
  }
  config[is.na(config)]="none"
  if(!dir.exists(tools::R_user_dir("iCRISPR")))dir.create(tools::R_user_dir("iCRISPR"),recursive = TRUE)
  config=data.frame(item=names(config),value=config)
  write.csv(config,file.path(tools::R_user_dir("iCRISPR"),"config"),quote = FALSE,row.names = F)
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
if(!dir.exists(tools::R_user_dir("iCRISPR")))dir.create(tools::R_user_dir("iCRISPR"),recursive = TRUE)
config=refresh_config()
