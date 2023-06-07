#' Summary Cas-system type
#'
#' @param crispr crispr object
#'
#' @export
#'
#' @examples
#' data(multi_crispr)
#' cas_type_res=summary_cas_type(multi_crispr)
#' plot(cas_type_res)
summary_cas_type=function(crispr){
  if(inherits(crispr,"multi_crispr")){
    aaa=lapply(crispr, summary_cas_type)
    bbb=do.call(rbind,aaa)
    if(is.null(bbb))return(NULL)
    rownames(bbb)=NULL
    class(bbb)=c("cas_type",class(bbb))
    return(bbb)
  }
  if(!inherits(crispr,"crispr"))return(NULL)
  cas_info=crispr$Cas
  if(is.null(cas_info))return(NULL)
  cas_type=cas_info[,c("genome","Cas_id","type","subtype")]%>%dplyr::distinct()%>%
    dplyr::count(genome,type,subtype)%>%as.data.frame()
  class(cas_type)=c("cas_type",class(cas_type))
  return(cas_type)
}

#' Plot the Cas system type and subtype
#' @method plot cas_type
#' @param x cas_type object
#' @param mode 1~3
#' @param ... additional
#'
#' @return ggplot
#' @exportS3Method
#' @rdname summary_cas_type
plot.cas_type=function(x,mode=1,...){
  cas_plotdat=dplyr::group_by(x[,c("type","subtype","n")],type,subtype)%>%dplyr::summarise(n=sum(n))%>%as.data.frame()
  if(mode==1)p=pcutils::gghuan2(cas_plotdat,...)
  if(mode==2)p=pcutils::gghuan2(cas_plotdat[,c("subtype","type","n")],...)
  if(mode==3)p=pcutils::my_sankey(cas_plotdat,mode = "gg",num=T,...)
  #if(mode==4)p=pcutils::my_circo(cas_plotdat,...)
  return(p)
}


#' Show cas type
#'
#' @return ggplot
#' @export
#'
#' @examples
#' show_cas_type()
#' @references 1. Makarova, K. S. et al. Evolutionary classification of CRISPR-Cas systems: a burst of class 2 and derived variants. Nat Rev Microbiol 18, 67–83 (2020).
show_cas_type=function(){
  pcutils::lib_ps("ggtree","gggenes","aplot",library = F)
  #load("cas_info.rda",envir = environment())
  data("cas_info.rda",package = "iCRISPR",envir = environment())
  #cas_type_df$cas1=strsplit(cas_type_df$cas,"[,-]")%>%sapply(.,\(i)i[1])

  cas_tree=ggtree::fortify(cas_tree)
  tree=ggtree::ggtree(cas_tree,ladderize = F)+
    ggtree::geom_tiplab()+
    ggtree::geom_highlight(node = 24,fill="skyblue",alpha=0.5,to.bottom = T)+
    ggtree::geom_highlight(node = 25,fill="pink",alpha=0.5,to.bottom = T)+
    geom_text(data = data.frame(x = c(1, 1),y = c(21.5, 9.5),
                                label = c("Class 1", "Class 2")),
              mapping = aes(x = x, y = y, label = label))+xlim(0,5)

  p=ggplot(cas_type_df,aes(y=sub_type,fill=cas1))+
    gggenes::geom_gene_arrow(aes(xmin=position,xmax=position+1,forward=(strand=="+")))+
    gggenes::geom_gene_label(aes(xmin=position,xmax=position+1,label=cas))+
    geom_hline(yintercept = 10.5,linetype=2)+
    theme_void()+theme(legend.position = "none")
  aplot::insert_left(p,tree,width = 0.3)
  #ggsave("inst/extdata/show_cas_type.pdf",width = 14,height = 6)
}

#' Summary crispr and spacer number at different evidence levels
#'
#' @param crispr crispr object
#' @param each_genome summary each genome information in multi_crispr?
#'
#' @export
#'
#' @examples
#' level_res=summary_levels(multi_crispr,each_genome=FALSE)
#' plot(level_res)
summary_levels=function(crispr,each_genome=T){
  if(inherits(crispr,"multi_crispr")){
    aaa=lapply(crispr, summary_levels)
    bbb=do.call(rbind,aaa)
    rownames(bbb)=NULL
    if(!each_genome){
      ccc=dplyr::group_by(bbb,evidence_level,Cas)%>%
        dplyr::summarise(crispr_num=sum(crispr_num),spacer_num=sum(spacer_num))%>%
        as.data.frame()
      ccc=dplyr::mutate(ccc,spacer_num_per_array=spacer_num/crispr_num)
      attributes(ccc)$each_genome=F
    }
    class(ccc)=c("evidence_level",class(ccc))
    return(ccc)
  }
  if(!inherits(crispr,"crispr"))return(NULL)
  if(is.null(crispr$Spacer))return(NULL)

  if("evidence_Level"%in%colnames(crispr$CRISPR))crispr$CRISPR=dplyr::rename(crispr$CRISPR,evidence_level="evidence_Level")

  if("spacer_number"%in%colnames(crispr$CRISPR)){
    n_spacer2=dplyr::select(crispr$CRISPR,genome,evidence_level,spacer_number)%>%
      dplyr::group_by(genome,evidence_level)%>%dplyr::summarise(crispr_num=dplyr::n(),spacer_num=sum(spacer_number))%>%
      as.data.frame()
  }
  else{
    #旧版本不包含spacer_number，需要统计
    n_spacer=dplyr::count(crispr$Spacer,genome,CRISPR_id)
    n_spacer2=dplyr::select(crispr$CRISPR,genome,CRISPR_id,evidence_level)%>%
      dplyr::left_join(.,n_spacer,by =c("genome"="genome","CRISPR_id"="CRISPR_id"))%>%
      dplyr::group_by(genome,evidence_level)%>%dplyr::summarise(crispr_num=dplyr::n(),spacer_num=sum(n))%>%as.data.frame()
  }
  n_spacer2$Cas=ifelse(is.null(crispr$Cas),"Non-Cas","With-Cas")

  class(n_spacer2)=c("evidence_level",class(n_spacer2))
  return(n_spacer2)
}


#' Plot the distribution of array and spacer number at different evidence levels.
#' @method plot evidence_level
#' @param x evidence_level object
#' @param mode 1~2
#' @param ... additional
#'
#' @return ggplot
#' @exportS3Method
#' @rdname summary_levels
plot.evidence_level=function(x,mode=1,...){
  pcutils::lib_ps("reshape2","scales",library = F)
  ccc=x
  if(attributes(x)$each_genome){
    ccc=dplyr::group_by(x,evidence_level,Cas)%>%
      dplyr::summarise(crispr_num=sum(crispr_num),spacer_num=sum(spacer_num))%>%
      as.data.frame()
    ccc=dplyr::mutate(ccc,spacer_num_per_array=spacer_num/crispr_num)
  }
  ccc$spacer_num_per_array=round(ccc$spacer_num_per_array,1)
  plotdat=reshape2::melt(ccc,id.vars = 1:2)
  if(mode==1){
    p=ggplot(data = plotdat,aes(x=evidence_level,y=value,fill=Cas))+
      geom_col(position = position_dodge(width = 1))+
      geom_text(aes(label=value),position = position_dodge(width = 1),vjust=0)+
      scale_fill_manual(values = c( "#a6cee3","#78c679"))+labs(y="Number")+
      facet_wrap(.~variable,nrow = 1,scales = "free_y")
  }
  if(mode==2){
    p=ggplot(data = plotdat,aes(x=evidence_level,y=value,fill=Cas))+
      geom_bar(stat="identity",position = position_fill())+
      scale_y_continuous(labels = scales::percent)+
      geom_text(aes(label=value),position = position_fill())+labs(y="Percentage")+
      scale_fill_manual(values = c( "#a6cee3","#78c679"))+
      facet_wrap(.~variable,nrow = 1)
  }
  p
}
