test_for_combine=function(func,all,parts){
  all_res=do.call(func,list(all))
  part_res=lapply(parts,func)
  part_res2=do.call(paste0("combine.",class(all_res)[1]),part_res)
  identical(all_res,part_res2)
}

#analysis after getting blast result

le=c("root","Kingdom","Phylum","Class","Order","Family","Genus","Species","Genome")

#' Summary the phylogenetic distance and hit count between source and target
#'
#' @param source_target_info data.frame which columns: "source_name","source_lineage","target_name","target_lineage"
#'
#' @return s_t_count
#' @export
#'
#' @examples
#' data(pro_net)
#' source_target_count(pro_net)
source_target_count=function(source_target_info){
  if(!all(c("source_name","source_lineage","target_name","target_lineage")%in%colnames(source_target_info)))stop("check columns!")
  arc_net=source_target_info
  #==count table
  arc_net%>%dplyr::select(source_name,target_name)%>%
    dplyr::count(source_name,target_name,name = "spacer_n")%>%data.frame()->s_t_count
  class(s_t_count)=c("s_t_count",class(s_t_count))
  s_t_count
}

#' Combine some s_t_count object
#'
#' @param ... s_t_count object
#'
#' @return s_t_count object
#' @exportS3Method
#' @method combine s_t_count
#'
combine.s_t_count=function(...){
  all_c=list(...)
  if(!all(sapply(all_c, inherits,what="s_t_count")))stop("all input should be s_t_count object")
  all_res=do.call(rbind,all_c)
  ccc=dplyr::group_by(all_res,source_name,target_name)%>%
    dplyr::summarise(spacer_n=sum(spacer_n))%>%
    as.data.frame()
  class(ccc)=c("s_t_count","data.frame")
  ccc
}

#' Summary the phylogenetic distance between source and target
#'
#' @param source_target_info data.frame which columns: "source_name","source_lineage","target_name","target_lineage"
#'
#' @return s_t_pd
#' @export
#'
#' @examples
#' data(pro_net)
#' source_target_pd(pro_net)
source_target_pd=function(source_target_info){
  if(!all(c("source_name","source_lineage","target_name","target_lineage")%in%colnames(source_target_info)))stop("check columns!")
  arc_net=source_target_info

  arc_net%>%dplyr::select(source_name,source_lineage,target_name,target_lineage)%>%
    dplyr::distinct_all()->s_t_pd

  pcutils::strsplit2(s_t_pd$source_lineage,";")%>%
    as.data.frame()->arc_taxonomy1
  pcutils::strsplit2(s_t_pd$target_lineage,";")%>%
    as.data.frame()->arc_taxonomy2
  res=data.frame(arc_taxonomy1==arc_taxonomy2)
  res1=apply(res,1,which)
  res2=sapply(res1, \(i)ifelse(length(i)==0,0,max(i)))
  s_t_pd=data.frame(s_t_pd[,c(1,3)],phy_dis=(7-res2)*2)
  class(s_t_pd)=c("s_t_pd",class(s_t_pd))
  s_t_pd
}

#' Combine some s_t_pd object
#'
#' @param ... s_t_pd object
#'
#' @return s_t_pd object
#' @exportS3Method
#' @method combine s_t_pd
#'
combine.s_t_pd=function(...){
  all_c=list(...)
  if(!all(sapply(all_c, inherits,what="s_t_pd")))stop("all input should be s_t_pd object")
  all_res=do.call(rbind,all_c)
  ccc=dplyr::distinct_all(all_res)
  class(ccc)=c("s_t_pd","data.frame")
  ccc
}

# 建树再计算比较慢
# source_target_pd=function(source_target_info, file = NULL){
#   if(!all(c("source_name","source_lineage","target_name","target_lineage")%in%colnames(source_target_info)))stop("check columns!")
#   arc_net=source_target_info
#   #==all species pd=
#   taxs=arc_net%>%dplyr::select(source_name,source_lineage)%>%dplyr::rename("target_name"=1,"target_lineage"=2)%>%
#     rbind(.,arc_net%>%dplyr::select(target_name,target_lineage))
#   taxs=taxs%>%dplyr::distinct(target_name,.keep_all = T)
#   pcutils::strsplit2(taxs$target_lineage,";",
#                      colnames =le[2:8] )%>%
#     as.data.frame()->arc_taxonomy
#   rownames(arc_taxonomy)<-taxs$target_name
#   pctax::df2tree(arc_taxonomy)->arc_tree
#   #treeio::as_tibble(arc_tree)->arc_tree2
#   stats::cophenetic(arc_tree)->p_dis
#
#   arc_net%>%dplyr::select(source_name,target_name)%>%
#     dplyr::distinct_all()->s_t_pd
#
#   apply(s_t_pd,1,\(x)p_dis[x[1],x[2]])->s_t_pd$phy_dis
#   if (!is.null(file)){
#     saveRDS(s_t_pd, file = file)
#     message(paste0("Result saved as ", file))
#   }
#   s_t_pd
# }


#' Summary the phylogenetic distance and hit count between source and target
#'
#' @param source_target_info data.frame which columns: "source_name","source_lineage","target_name","target_lineage"
#'
#' @return s_t_gng
#' @export
#'
#' @examples
#' data(pro_net)
#' source_target_gng(pro_net)
source_target_gng=function(source_target_info){
  if(!all(c("source_name","source_lineage","target_name","target_lineage","gene_or_not")%in%colnames(source_target_info)))stop("check columns!")
  arc_net=source_target_info
  #==count table
  gene_count=dplyr::count(arc_net,source_name,target_name,gene_or_not)%>%
    reshape2::dcast(source_name+target_name~gene_or_not,value.var="n")
  gene_count[is.na(gene_count)]=0
  s_t_gng=dplyr::mutate(gene_count,gene_ratio=gene/(gene+intergene))
  class(s_t_gng)=c("s_t_gng",class(s_t_gng))
  s_t_gng
}

#' Combine some s_t_gng object
#'
#' @param ... s_t_gng object
#'
#' @return s_t_gng object
#' @exportS3Method
#' @method combine s_t_gng
#'
combine.s_t_gng=function(...){
  all_c=list(...)
  if(!all(sapply(all_c, inherits,what="s_t_gng")))stop("all input should be s_t_gng object")
  all_res=do.call(rbind,all_c)
  ccc=dplyr::group_by(all_res,source_name,target_name)%>%
    dplyr::summarise(gene=sum(gene),intergene=sum(intergene))%>%
    as.data.frame()
  ccc=dplyr::mutate(ccc,gene_ratio=gene/(gene+intergene))
  class(ccc)=c("s_t_gng","data.frame")
  ccc
}

gene_pd_lm=function(arc_count2,mode=1,filter_spacer=5){
  arc_count2=dplyr::filter(arc_count,spacer_n>filter_spacer)
  if(mode==1){
    p<-ggplot(arc_count2,aes(phy_dis,gene_ratio))
    p <-p+geom_jitter(size = 0.1,alpha=0.5) +
      geom_smooth(method = "lm",
                  color = "red", se = F, formula = "y~x") +
      ggpmisc::stat_poly_eq(aes(label = paste(
        #after_stat(eq.label),
        after_stat(adj.rr.label), after_stat(p.value.label),sep = "~~~~~")),
        formula = y ~ x, parse = TRUE, color = "blue", size = 3, label.x = 0.05, label.y = 1.2) + labs(x = NULL, y = NULL)+
      #stat_cor(method = "pearson",color='red')+
      theme_classic() +
      scale_y_continuous(breaks = seq(0,1,0.2))+
      labs(x='Phylogenetic distance',y='coding proportion') +
      ggpubr::theme_pubr(base_size=20)+
      theme(axis.text = element_text(size=16))
  }
  if(mode==2){
    p <-ggplot(arc_count2,aes(phy_dis,gene_ratio))+
      geom_boxplot(aes(phy_dis,gene_ratio,group=phy_dis),outlier.shape = NA)+
      geom_jitter(size = 0.3,alpha=0.8)+
      #geom_hline(yintercept = mean_value, linetype = "dashed", color = "red")+
      geom_smooth(method = "lm",
                  color = "red", se = F, formula = "y~x") +
      ggpmisc::stat_poly_eq(aes(label = paste(
        #after_stat(eq.label),
        after_stat(adj.rr.label), after_stat(p.value.label),sep = "~~~~~")), formula = y ~
          x, parse = TRUE, color = "blue", size = 3, label.x = 0.05,
        label.y = 1.2) + labs(x = NULL, y = NULL)+
      scale_x_continuous(breaks=seq(0, 14, 2))+
      scale_y_continuous(breaks = seq(0,1,0.2))+
      labs(x='Phylogenetic distance(Archaea)',y='coding proportion')+
      ggpubr::theme_pubr(base_size=20)+
      theme(axis.text = element_text(size=16))
  }
  p
}


#' Sankey plot for source-target
#'
#' @param source_target_info data.frame which columns: "source_lineage","target_lineage"
#' @param two_level manual which two levels in "Kingdom","Phylum","Class","Order","Family","Genus","Species". Sometimes "Genome"
#' @param topN1 select topN for level1
#' @param topN2 select topN for level2
#' @param file filename
#'
#' @return ggplot
#' @export
#'
#' @examples
#' sankey_overview(pro_net,two_level=c("Species","Genus"))
sankey_overview=function(source_target_info,two_level=NULL,topN1=8,topN2=8,file=NULL){
  if(!all(c("source_lineage","target_lineage")%in%colnames(source_target_info)))stop("check columns!")
  arc_net=source_target_info
  pcutils::strsplit2(arc_net$source_lineage,";")%>%as.data.frame()->source_tax
  pcutils::strsplit2(arc_net$target_lineage,";")%>%as.data.frame()->target_tax
  colnames(source_tax)=paste0("source_",le[2:8])
  if("source_genome"%in%colnames(arc_net))source_tax$source_Genome=arc_net$source_genome
  colnames(target_tax)=paste0("target_",le[2:8])
  if("target_genome"%in%colnames(arc_net))source_tax$target_Genome=arc_net$target_genome
  cbind(source_tax,target_tax)->source_target

  if(!is.null(two_level)){
    two_level=rep(two_level,len=2)
    source_target%>%dplyr::select(paste0("source_",two_level[1]),paste0("target_",two_level[2]))%>%
      MetaNet::summ_2col(direct = T)%>%dplyr::arrange(-count)->a
    p=pcutils::my_sankey(a,mode = "gg",topN = topN1,num=T)
    return(p)
  }
  #前10的link
  pls=list()
  for (levels in le[2:8]) {
    source_target%>%dplyr::select(dplyr::ends_with(levels))%>%MetaNet::summ_2col(direct = T)%>%dplyr::arrange(-count)->a
    pls[[levels]]=pcutils::my_sankey(a,mode = "gg",topN = topN1,num=T)
    #看一个source tax的所有target
    pls1=list()
    for (i in unique(a[,1])[1:10]) {
      a2=a[a[,1]==i,]%>%as.data.frame()
      pls1[[i]]=pcutils::my_sankey(a2,mode = "gg",topN = topN2,num=F)
    }
    pcutils::plotpdf(pls1,paste0(file,"_",levels,"_single_source"),height = 10)
  }
  pcutils::plotpdf(pls,paste0(file,"each_level"),height = 10)
  message("All done.")
}

#' Summary spacers target which: Self? Virus? Others?
#'
#' @param source_target_info data.frame which columns: "source_lineage","target_lineage"
#' @param two_level manual which two levels in "Kingdom","Phylum","Class","Order","Family","Genus","Species". Sometimes "Genome"
#'
#' @return data.frame
#' @export
#'
#' @examples
#' tri_target(pro_net,two_level=c("Species","Genus"))
tri_target=function(source_target_info,two_level=NULL){
  if(!all(c("source_lineage","target_lineage")%in%colnames(source_target_info)))stop("check columns!")
  arc_net=source_target_info
  pcutils::strsplit2(arc_net$source_lineage,";")%>%as.data.frame()->source_tax
  pcutils::strsplit2(arc_net$target_lineage,";")%>%as.data.frame()->target_tax
  colnames(source_tax)=paste0("source_",le[2:8])
  if("source_genome"%in%colnames(arc_net))source_tax$source_Genome=arc_net$source_genome
  colnames(target_tax)=paste0("target_",le[2:8])
  if("target_genome"%in%colnames(arc_net))source_tax$target_Genome=arc_net$target_genome
  cbind(source_tax,target_tax)->source_target

  get_self_target=function(source_target,levels){
    if(which(levels[1]==le)<which(levels[2]==le))stop("levels[1] should higher than levels[2] likes c(\"Species\",\"Genus\")")
    dat=data.frame("id"=dplyr::pull(source_target,paste0("source_",levels[1])),
                 "source"=dplyr::pull(source_target,paste0("source_",levels[2])),
                 "target"=dplyr::pull(source_target,paste0("target_",levels[2])),
                 "kingdom"=dplyr::pull(source_target,target_Kingdom))

    a=dat%>%dplyr::group_by_all()%>%dplyr::count()%>%as.data.frame()
    #self-targeting fraction
    tmpls=lapply(unique(a[,"id"]),\(i){
      tmpdf=a[a[,"id"]==i,]%>%as.data.frame()
      source=unique(tmpdf[,"source"])
      if(length(source)>1)stop(i," has multiple ",levels[2]," !")
      all=sum(tmpdf[,"n"])
      s=tmpdf[tmpdf[,"target"]==source,"n"]
      if(length(s)==0)s=0
      else s=sum(s)
      v=tmpdf[tmpdf[,"kingdom"]=="k__Viruses","n"]
      if(length(v)==0)v=0
      else v=sum(v)
      o=all-s-v
      return(data.frame("Source"=i,"Spacer_count"=all,"T_self"=s,"T_virus"=v,"T_others"=o))
    })
    self_df=do.call(rbind,tmpls)
    self_df$`T_self_r`=round(self_df$`T_self`/self_df$`Spacer_count`,4)
    self_df$`T_virus_r`=round(self_df$`T_virus`/self_df$`Spacer_count`,4)
    self_df$`T_others_r`=round(1-self_df$T_self_r-self_df$T_virus_r,4)
    if(levels[1]==levels[2])levels=levels[1]
    else levels=paste(levels,collapse = "-")
    data.frame(Level=levels,self_df)
  }

  if(!is.null(two_level)){
    two_level=rep(two_level,len=2)
    return(get_self_target(source_target,two_level))
  }
  res=lapply(le[2:8], \(i){
    levels=rep(i,len=2)
    get_self_target(source_target,levels)
  })
  self_df=do.call(rbind,res)
  self_df$Level=pcutils::change_fac_lev(self_df$Level,le[2:8])
  rownames(self_df)=self_df$Source
  class(self_df)=c("self_df",class(self_df))
  self_df
}

#' Combine some self_df object
#'
#' @param ... self_df object
#'
#' @return self_df object
#' @exportS3Method
#' @method combine self_df
#'
combine.self_df=function(...){
  all_c=list(...)
  if(!all(sapply(all_c, inherits,what="self_df")))stop("all input should be self_df object")
  all_res=do.call(rbind,all_c)
  self_df=dplyr::group_by(all_res,Level,Source)%>%
    dplyr::summarise(Spacer_count=sum(Spacer_count),T_self=sum(T_self),
                     T_virus=sum(T_virus),T_others=sum(T_others))%>%
    as.data.frame()

  self_df$`T_self_r`=round(self_df$`T_self`/self_df$`Spacer_count`,4)
  self_df$`T_virus_r`=round(self_df$`T_virus`/self_df$`Spacer_count`,4)
  self_df$`T_others_r`=round(1-self_df$T_self_r-self_df$T_virus_r,4)
  self_df$Level=pcutils::change_fac_lev(self_df$Level,le[2:8])
  self_df=dplyr::arrange(self_df,Level,Source)
  rownames(self_df)=self_df$Source
  class(self_df)=c("self_df","data.frame")
  self_df
}

#' plot self_df
#'
#' @param x self_df
#' @param mode 1~2
#' @param filter_spacer filter spacer less than filter_spacer, default 0
#' @param ... add
#'
#' @return ggplot
#' @exportS3Method
#' @method plot self_df
plot.self_df=function(x,mode=1,filter_spacer=0,...){
  self_df=x
  self_df%>%filter(Spacer_count>filter_spacer,Level!="Kingdom")->tmp
  if(mode==1){
    p = ggtern::ggtern(tmp, aes(x = T_self_r, y = T_virus_r, z = T_others_r)) +
      geom_point(aes(size = Spacer_count, col = Level)) +ggsci::scale_color_d3()+
      theme_bw()+  ggtern::theme_showarrows()+labs(title = paste0("Fitler spacer number >",filter_spacer))
    }

  if(mode==2){
    p = ggtern::ggtern(tmp, aes(x = T_self_r, y = T_virus_r, z = T_others_r)) +
      geom_point(aes(size = Spacer_count, col = Level)) +theme_bw()+ ggtern::theme_showarrows()+
      ggsci::scale_color_d3()+facet_wrap(.~Level)+labs(title = paste0("Fitler spacer number >",filter_spacer))
    }
  return(p)
}


#' Build tri_target_tree
#'
#' @param source_target_info data.frame which columns: "source_lineage","target_lineage"
#' @param self_df self_df from `tri_target`
#' @param level 1~7, means "Kingdom","Phylum","Class","Order","Family","Genus","Species"
#' @param filter_spacer filter spacer less than filter_spacer, default 0
#'
#' @return tree
#' @export
#'
#' @examples
#' self_df=tri_target(pro_net)
#' tree=tri_target_tree(pro_net,self_df,filter_spacer=10)
#' #plot(tree)
tri_target_tree=function(source_target_info,self_df,level=7,filter_spacer=0){
  if(!all(c("source_lineage","target_lineage")%in%colnames(source_target_info)))stop("check columns!")
  arc_net=source_target_info
  pcutils::strsplit2(arc_net$source_lineage,";",colnames =le[2:8])%>%
    as.data.frame()->source_tax
  dplyr::distinct_all(source_tax)->source_lineage

  pcutils::lib_ps("ggtree","treedataverse",library = FALSE)

  f_tax=distinct(source_lineage[,1:level])
  rownames(f_tax)=f_tax[,level]
  self_df_f=filter(self_df,Spacer_count>filter_spacer,Level==le[level+1])
  f_tax=f_tax[rownames(self_df_f),]
  otutab=self_df[rownames(self_df_f),3,drop=F]

  pctax::ann_tree(f_tax,otutab,level = level)->tree

  tree=left_join(tree,self_df,by=c("label"="Source"))

  tree%>%mutate(label1=sub(".?__","",label))->tree2
  class(tree2)=c("tri_target_tree",class(tree2))
  attributes(tree2)$filter_spacer=filter_spacer
  tree2
}

#' Plot tri_target_tree
#'
#' @param x tri_target_tree
#' @param add_self logical
#' @param add_virus logical
#' @param add_others logical
#' @param add_tiplab logical
#' @param some_tax vectors
#' @param ... add
#' @param label_size label size
#'
#' @return ggplot
#' @exportS3Method
#' @method plot tri_target_tree
plot.tri_target_tree=function(x,add_self=TRUE,add_virus=TRUE,add_others=TRUE,add_tiplab=TRUE,
                              some_tax=c(),label_size=1.5,...){
  tree=x
  p=pctax::easy_tree(tree,add_abundance = F,add_tiplab = F,...)
  if(add_self){
    p=p+ggnewscale::new_scale_fill()+
      ggtreeExtra::geom_fruit(
        geom=geom_tile,
        mapping=aes(y=label, fill=T_self_r),
        stat="identity",offset = 0.12,pwidth = 0.1,
      )+scale_fill_gradient(low = "#FFFFFF",high = "blue")
  }
  if(add_virus){
    p=p+ggnewscale::new_scale_fill()+
      #scale_fill_d3()+
      ggtreeExtra::geom_fruit(
        geom=geom_tile,
        mapping=aes(y=label, fill=T_virus_r),
        stat="identity",offset = 0.12,pwidth = 0.1,
      )+scale_fill_gradient(low = "#FFFFFF",high = "red")
  }
  if(add_virus){
    p=p+ggnewscale::new_scale_fill()+
      ggtreeExtra::geom_fruit(
        geom=geom_tile,
        mapping=aes(y=label, fill=T_others_r),
        stat="identity",offset = 0.12,pwidth = 0.1,
    )+scale_fill_gradient(low = "#FFFFFF",high = "green")
  }
  if(add_tiplab){
    p=p+ggtree::geom_tiplab(aes(label=label1),color="black",size=label_size,
                            offset = 3, show.legend = FALSE)
  }
  if(length(some_tax)>0){
    p=pctax::add_strip(p,some_tax = some_tax,strip_params = list(offset=1))
  }
  p+ggtitle(paste0("Fitler spacer number >",attributes(tree)$filter_spacer))
}

# data2<-filter(arc_net2, target_kingdom%in%c("Archaea","Viruses"))
#
# p <- ggplot(data2, aes(relative_loci,..density..,fill=target_kingdom,color=target_kingdom))+
#   geom_histogram(binwidth = 0.01,alpha=0.5)+
#   geom_line(stat='density',size=0.8,color="black")+
#   theme_bw() + theme(panel.grid=element_blank())+scale_y_continuous(expand=c(0,0))+
#   scale_color_manual(values = c("#E9967A","#483D8B"))+
#   scale_fill_manual(values = c("#E9967A","#483D8B"))+#"#838B8B" cellular organism;"#B2DFEE" protozoa
#   #ggtitle("relative position in gene")+
#   xlab("Protospacers' relative position in genes")+
#   scale_x_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1))+
#   ggpubr::theme_pubr(base_size = 24,legend = "none")+facet_wrap(target_kingdom~.)
#
# p
#
# ggsave(filename = "analysis/arc23/arc_targets_relative_pos.pdf",
#        plot=p,
#        device = "pdf",
#        width = 8,
#        height = 4,
#        dpi = 600)
