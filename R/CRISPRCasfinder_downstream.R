#' Prepare the result files from CRISPRCasFinder
#'
#' @param input_folder the folder of CRISPRCasFinder result, contains GFF/,TSV/,result.json,...
#' @param output_folder output, default: ./
#' @param genome_name the genome_name
#'
#' @export
#'
#' @examples
#' crispr=pre_CCF_res(system.file("extdata/AAB-S01R1_122",package = "iCRISPR"))
pre_CCF_res=function(input_folder,output_folder="./pre_CCF_res_out",genome_name=NULL){
  if(!dir.exists(input_folder))stop("Can not find ",input_folder)
  if(!dir.exists(output_folder))dir.create(output_folder,recursive =T)

  if(is.null(genome_name))genome_name=basename(input_folder)

  out_file1=paste0(output_folder,"/",genome_name,"_spacer.fa")
  out_file2=paste0(output_folder,"/",genome_name,"_spacer_with_cas.fa")

  #有annotation前缀的gff文件是有crispr array的contigs
  crispr_gff_ls=list.files(paste0(input_folder,"/GFF"),pattern = "annotation*")%>%
    gsub("annotation_","",.)

  if(length(crispr_gff_ls)>0){
    arrary_res=res=data.frame()
    for (i in crispr_gff_ls) {
      crispr_gff=paste0(input_folder,"/GFF/",i)
      tmp=get_spacer_fa(crispr_gff,genome_name=genome_name)
      res=rbind(res,tmp)
      tmp=get_array(crispr_gff,genome_name=genome_name)
      arrary_res=rbind(arrary_res,tmp)
    }
    write_fasta(res,out_file1)

    #导出Spacer信息
    write.csv(arrary_res,file = paste0(output_folder,"/",genome_name,"_Array_info.csv"),row.names = F)
    spacer_res=filter(arrary_res,feature=="CRISPRspacer")
    write.csv(spacer_res,file = paste0(output_folder,"/",genome_name,"_Spacer_info.csv"),row.names = F)

    crispr_res=get_crispr(paste0(input_folder,"/TSV/Crisprs_REPORT.tsv"),genome_name = genome_name)
    #导出CRISPR信息
    write.csv(crispr_res,file = paste0(output_folder,"/",genome_name,"_CRISPR_info.csv"),row.names = F)
  } else {
    pcutils::dabiao("No array found!")
  }

  #导出Cas信息
  rawCas_file=paste0(input_folder,"/rawCas.fna")
  Cas_report_file=paste0(input_folder,"/TSV/Cas_REPORT.tsv")
  if(file.exists(rawCas_file)){
    cas_res=get_cas(rawCas_file,Cas_report_file,genome_name = genome_name)
    write.csv(cas_res,file = paste0(output_folder,"/",genome_name,"_Cas_info.csv"),row.names = F)
  }else {
    cas_res=NULL
    pcutils::dabiao("No Cas found!")
  }

  #导出包含cas的spacer
  # df=read_fasta(out_file1)
  # spacer_with_cas=df[grepl(paste0(cas_res$seqid,collapse = "|"),df$Sequence_ID),]
  # write_fasta(spacer_with_cas,out_file2)

  pcutils::dabiao("All done! see ",output_folder,"/")
  return(list(CRISPR=crispr_res,Cas=cas_res,Arrary=arrary_res,Spacer=spacer_res))
}

#可能不太正确，因为一个genome是环状的，存在Cas有可能全部都有
get_spacer_fa=function(crispr_gff,genome_name){
  #用@分割，header要包括这些信息：spacer来源于哪条genome，哪条序列，是在该序列第几个crispr array 上，crispr array的起始和终止位点，再加一个spacer的起始位点，spacer的长度，以及是第几个spacer
  res=data.frame(spacer_id=c(),spacer_seq=c())
  gff1=read.table(crispr_gff,col.names = c("seqid", "source", "feature", "start", "end", "score", "strand", "phase", "attributes"))
  gff2=filter(gff1,feature%in%c("CRISPR","CRISPRspacer"))%>%select(1,3,4,5,9)
  cid=0
  for (r in 1:nrow(gff2)) {
    row=gff2[r,]
    if(row$feature=="CRISPR"){
      cid=cid+1
      sid=0
      crispr1=row
    }
    else if(row$feature=="CRISPRspacer") {
      sid=sid+1
      attri=strsplit(row$attributes,";")[[1]]
      spacer_id=paste0(genome_name,
                       "@",row$seqid,
                       "@CRISPR:",cid,
                       "@",crispr1$start,"-",crispr1$end,
                       "@",gsub("Name=","",attri[2]),"_",sid)
      spacer_seq=gsub("sequence=","",attri[1])
      # if(is.null(res))res=paste0(spacer_id,"\n",spacer_seq)
      # else res=paste0(res,"\n",spacer_id,"\n",spacer_seq)
      res=rbind(res,data.frame(spacer_id=spacer_id,spacer_seq=spacer_seq))
    }
  }
  return(res)
}

get_array=function(crispr_gff,genome_name){
  gff=read.file(crispr_gff)
  left=which(gff$feature=="LeftFLANK")
  right=which(gff$feature=="RightFLANK")
  if(length(left)!=length(right))stop("Something wrong with this gff file: ",crispr_gff)
  crispr_id=c()
  for (i in seq_along(left)) {
    crispr_id[left[i]:right[i]]=paste0("CRISPR:",i)
  }
  gff$CRISPR_id=crispr_id
  gff$CRISPR_id=paste0(gff$seqid,"@",gff$CRISPR_id)
  gff=filter(gff,feature!="CRISPR")
  array_res=gff[,c(1,10,3:5,7)]
  array_res$sequence=sub(".*sequence=([^;]+).*","\\1",gff$attributes)
  array_res$strand=ifelse(array_res$strand=="+","Forward",ifelse(array_res$strand=="-","Reverse","Unknown"))
  return(data.frame(genome=genome_name,array_res))
}

# get_cas=function(rawCas_file, Cas_report_file,genome_name){
#   rawCas_fasta=read_fasta(rawCas_file)%>%distinct(Sequence_ID,.keep_all = T)
#   rawCas_fasta$Sequence_ID=lapply(rawCas_fasta$Sequence_ID,\(i)strsplit(i,split="\\|")[[1]][2])%>%do.call("c",.)
#   rawCas_fasta=rawCas_fasta%>%distinct(Sequence_ID,.keep_all = T)
#
#   cas_report=readLines(Cas_report_file)
#   cas_report=cas_report[!(grepl("^[#-]",cas_report)|grepl("^$",cas_report)|grepl("\\( Unknown \\)",cas_report))]%>%
#     read.table(text = .,check.names = F,
#                col.names = c("Sequence_ID", "Cas-type/subtype", "Gene status", "System", "Type", "Begin", "End", "Strand", "Other_information")
# )
#   cas_report=cas_report%>%select("Sequence_ID","Type", "System","Cas-type/subtype", "Begin", "End", "Strand")%>%distinct()
#
#   cas_res=left_join(cas_report,rawCas_fasta,by=c("Sequence_ID"="Sequence_ID"))
#   cas_res$Type="Cas"
#   colnames(cas_res)=c("seqid","feature","subtype","protein","start","end","strand","sequence")
#   cas_res$strand=ifelse(cas_res$strand=="+","Forward",ifelse(cas_res$strand=="-","Reverse","Unknown"))
#   cas_res$seqid=gsub("_\\d+$","",cas_res$seqid)
#
#   cas_res$subtype=gsub("CAS-","",cas_res$subtype)
#   cas_res$type=ifelse(grepl("Type",cas_res$subtype),substr(cas_res$subtype,1,nchar(cas_res$subtype)-1),cas_res$subtype)
#
#   cas_res=cbind(genome=genome_name,cas_res[,c("seqid","feature" ,"type","subtype","protein","start","end","strand","sequence" )])
#   cas_res
# }

get_cas=function(rawCas_file, Cas_report_file,genome_name){
  all_cas_sys=data.frame()
  cas_report=readLines(Cas_report_file)
  cas_sys=c()
  flag=F
  flag2=F
  contig=""

  for (i in cas_report) {

    #System开始部分
    if(grepl("### System:",i)){
      sys_type=sub("### System: (.*) \\(.*\\)","\\1",i)
      flag=T
    }
    #System结束部分
    else if(grepl("####Summary system",i)){
      flag=F
      a=read.table(text = cas_sys,check.names = F,comment.char = "#",
                 col.names = c("Sequence_ID", "Cas-type/subtype", "Gene status", "System", "Type", "Begin", "End", "Strand", "Other_information")
      )
      a=a%>%select("Sequence_ID","Type", "System","Cas-type/subtype", "Begin", "End", "Strand")%>%distinct()

      contig1=sub(".*sequenceID=(\\S+)}.*","\\1",i)
      if(contig1!=contig){
        cas_id=1
      }
      else cas_id=cas_id+1
      contig=contig1

      a$Cas_id=paste0("CAS:",cas_id)
      a$Type=sys_type
      cas_sys=c()
      all_cas_sys=rbind(all_cas_sys,a)
    }
    #System开始时读入新行
    if(flag)cas_sys=c(cas_sys,i)
  }

  rawCas_fasta=read_fasta(rawCas_file)%>%distinct(Sequence_ID,.keep_all = T)
  rawCas_fasta$Sequence_ID=lapply(rawCas_fasta$Sequence_ID,\(i)strsplit(i,split="\\|")[[1]][2])%>%do.call("c",.)
  rawCas_fasta=rawCas_fasta%>%distinct(Sequence_ID,.keep_all = T)
  cas_res=left_join(all_cas_sys,rawCas_fasta,by=c("Sequence_ID"="Sequence_ID"))
  cas_res$feature="Cas"
  colnames(cas_res)=c("seqid","type","subtype","protein","start","end","strand","Cas_id","sequence","feature")
  cas_res$strand=ifelse(cas_res$strand=="+","Forward",ifelse(cas_res$strand=="-","Reverse","Unknown"))
  cas_res$seqid=gsub("_\\d+$","",cas_res$seqid)

  cas_res$subtype=gsub("CAS-","",cas_res$type)
  cas_res$type=ifelse(grepl("Type",cas_res$subtype),substr(cas_res$subtype,1,nchar(cas_res$subtype)-1),cas_res$subtype)

  cas_res$Cas_id=paste0(cas_res$seqid,"@",cas_res$Cas_id)
  cas_res=cbind(genome=genome_name,cas_res[,c("seqid","Cas_id","feature" ,"type","subtype","protein","start","end","strand","sequence" )])
  cas_res
}

get_crispr=function(Crisprs_REPORT,genome_name = genome_name){
  #Crisprs_REPORT="extdata/AAB-S01R1_122/TSV/Crisprs_REPORT.tsv"
  crisprs=read.table(Crisprs_REPORT,sep="\t",header = T,check.names = F)
  crisprs=crisprs[,c("Sequence","CRISPR_Id","Strain","CRISPR_Start","CRISPR_End","Potential_Orientation (AT%)","Consensus_Repeat")]

  crisprs$Strain="CRISPR"
  crisprs$CRISPR_Id=paste0(crisprs$Sequence,"@CRISPR:",sub(".*_(\\d+)$","\\1",crisprs$CRISPR_Id))
  colnames(crisprs)=c("seqid","CRISPR_id","feature","start","end","strand","consensus_repeat")
  data.frame(genome=genome_name,crisprs)
}

#' Plot a crispr-cas system
#'
#' @param crispr crispr result from `pre_CCF_res()`
#' @param genome the genome name
#' @param contig the contig name
#'
#' @return ggplot
#' @export
#'
#' @examples
#' plot_crispr(crispr,genome="AAB-S01R1_122",contig="AAB-S01R1_k55_9399631_flag=0_multi=9.8751_len=26518")
plot_crispr=function(crispr,genome=NULL,contig=NULL){
  #尝试画一下，类似基因图
  lib_ps("ggnewscale","gggenes","ggpubr")

  cas_res=crispr$Cas%>%filter(seqid==contig)
  cas_res=mutate(cas_res,protein=gsub("_.*","",protein))
  crispr_res=crispr$CRISPR%>%filter(seqid==contig)
  array_res=crispr$Arrary%>%filter(seqid==contig)
  array_res$feature=factor(array_res$feature,levels = c("LeftFLANK","CRISPRdr","CRISPRspacer","RightFLANK"))%>%droplevels()

  # ggplot()+
  #   geom_point(data = cas_res,aes(x=start,y=end),col="red")+
  #   geom_point(data = filter(crispr_res,seqid%in%cas_res$seqid),aes(x=start,y=end),col="blue")+
  #   scale_x_continuous(breaks = scales::pretty_breaks(n = 10))+
  #   facet_wrap(.~seqid,scales = "free",ncol = 1)

  p=ggplot()
  if(nrow(cas_res)>0){
    p=p+
      #cas基因
      geom_gene_arrow(data = cas_res, aes(xmin = start, xmax = end,y = seqid, fill = protein),
                      arrowhead_width=grid::unit(5,"mm"),arrowhead_height = grid::unit(5,"mm"),arrow_body_height= grid::unit(4,"mm")) +
      geom_gene_label(data = cas_res, aes(xmin = start, xmax = end,y = seqid, fill = protein,label=protein),min.size=0)+
      ggnewscale::new_scale_fill()
  }
  if(nrow(array_res)>0){
    p=p+
      #crispr系统
      geom_gene_arrow(data = array_res, aes(xmin = start, xmax = end,y = seqid,fill=feature),color=NA,
                      arrowhead_width=grid::unit(0,"mm"),arrowhead_height = grid::unit(0,"mm"))+
      scale_fill_manual(values = setNames(c("#5950FA","#190861","#8dd3c7","#fb8072"),
                                          c("LeftFLANK","CRISPRdr","CRISPRspacer","RightFLANK")))
  }
  p=p+
    coord_fixed(diff(range(c(cas_res$start,array_res$start,cas_res$end,array_res$end)))/8)+
    #ggpubr::theme_pubr()+
    theme_genes()+
    labs(y=NULL,title = paste0("Genome: ",genome," Contig: ",contig))+
    theme(legend.position = "top",
          axis.line.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
  p
}

#' Plot the Cas system type and subtype
#'
#' @param crispr crispr result from `pre_CCF_res()`
#'
#' @return ggplot
#' @export
#'
#' @examples
#' plot_cas_type(crispr)
plot_cas_type=function(crispr){
  cas_info=crispr$Cas
  a=cas_info[,c("Cas_id","type","subtype")]%>%distinct()%>%group_by_all()%>%
    count()%>%arrange(type)%>%as.data.frame()

  gghuan2(a[,c(2,3,4)])
}
