#' Prepare the result files from CRISPRCasFinder
#'
#' @param input_folder the folder of CRISPRCasFinder result, contains GFF/,TSV/,result.json,...
#' @param output_folder output, default: ./
#' @param genome_name the genome_name
#'
#' @export
#'
#' @examples
#' pre_CCF_res("extdata/AAB-S01R1_122","extdata")
pre_CCF_res=function(input_folder,output_folder="./",genome_name=NULL){
  if(!dir.exists(input_folder))stop("Can not find ",input_folder)
  if(!dir.exists(output_folder))dir.create(output_folder)

  if(is.null(genome_name))genome_name=basename(input_folder)

  out_file1=paste0(output_folder,"/",genome_name,"_spacer.fa")
  out_file2=paste0(output_folder,"/",genome_name,"_spacer_with_cas.fa")

  #有annotation前缀的gff文件是有crispr array的contigs
  crispr_gff_ls=list.files(paste0(input_folder,"/GFF"),pattern = "annotation*")%>%
    gsub("annotation_","",.)

  if(length(crispr_gff_ls)>0){
    res=data.frame()
    for (i in crispr_gff_ls) {
      crispr_gff=paste0(input_folder,"/GFF/",i)
      tmp=get_spacer(crispr_gff,genome_name=genome_name)
      res=rbind(res,tmp)
    }
    write_fasta(res,out_file1)

    crispr_res=get_crispr(rawCas_file,genome_name = genome_name)
    write.csv(crispr_res,file = paste0(output_folder,"/",genome_name,"_CRISPR_info.csv"),row.names = F)
  } else {
    pcutils::dabiao("No array found!")
  }

  #导出Cas信息
  rawCas_file=paste0(input_folder,"/rawCas.fna")
  if(file.exists(rawCas_file)){
    cas_res=get_cas(rawCas_file,genome_name = genome_name)
    write.csv(cas_res,file = paste0(output_folder,"/",genome_name,"_Cas_info.csv"),row.names = F)
  }else {
    pcutils::dabiao("No Cas found!")
  }

  #导出包含cas的spacer
  df=read_fasta(out_file1)
  spacer_with_cas=df[grepl(paste0(cas_res$seqid,collapse = "|"),df$Sequence_ID),]
  write_fasta(spacer_with_cas,out_file2)
  pcutils::dabiao("All done! see ",output_folder,"/")
}

get_spacer=function(crispr_gff,genome_name){
  #用@分割，header要包括这些信息：spacer来源于哪条genome，哪条序列，是在该序列第几个crispr array 上，crispr array的起始和终止位点，再加一个spacer的起始位点，spacer的长度，以及是第几个spacer
  res=data.frame(spacer_id=c(),spacer_seq=c())
  gff1=read.table(crispr_gff,col.names = c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"))
  gff2=filter(gff1,type%in%c("CRISPR","CRISPRspacer"))%>%select(1,3,4,5,9)
  cid=0
  for (r in 1:nrow(gff2)) {
    row=gff2[r,]
    if(row$type=="CRISPR"){
      cid=cid+1
      sid=0
      crispr1=row
    }
    else if(row$type=="CRISPRspacer") {
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

get_cas=function(rawCas_file,genome_name){
  #输出含有 cas protein 的 array
  arrary_with_cas_id=readLines(rawCas_file)%>%grep("^>",.,value = T)
  arrary_with_cas=pcutils::strsplit2(arrary_with_cas_id,"\\|")
  arrary_cas=pcutils::strsplit2(arrary_with_cas$V4,"[ ,]")
  res=cbind(arrary_with_cas[,c(2,3)],arrary_cas)
  colnames(res)=c("seqid","subtype","protein","start","end")
  res$seqid=gsub("_\\d+$","",res$seqid)
  res=dplyr::distinct(res)
  res$subtype=gsub("CAS-","",res$subtype)
  res$type=ifelse(grepl("Type",res$subtype),substr(res$subtype,1,nchar(res$subtype)-1),res$subtype)
  cas_res=cbind(genome=genome_name,res[,c("seqid","type","subtype","protein","start","end")])
  cas_res=mutate_at(cas_res,6:7,as.numeric)
  cas_res
}

get_crispr=function(Crisprs_REPORT,genome_name = genome_name){
  Crisprs_REPORT="extdata/AAB-S01R1_122/TSV/Crisprs_REPORT.tsv"
  crisprs=read.table(Crisprs_REPORT,sep="\t",header = T,check.names = F)
  crisprs=crisprs[,c("Sequence","CRISPR_Id","CRISPR_Start","CRISPR_End","Consensus_Repeat")]

  crisprs$CRISPR_Id=paste0(crisprs$Sequence,"@CRISPR:",sub(".*_(\\d+)$","\\1",crisprs$CRISPR_Id))
  colnames(crisprs)=c("seqid","CRISPR_id","start","end","consensus_repeat")
  data.frame(genome_name=genome_name,crisprs)
}


if(F){
  cas_info=read.csv("extdata/AAB-S01R1_122_Cas_info.csv")
  a=cas_info[,c("type","subtype")]%>%group_by_all()%>%count()%>%as.data.frame()

  my_sankey(a,mode = "gg")
  gghuan2(a[,c(1,2,3)])

  ggplot()+
    geom_point(data = cas_res,aes(x=start,y=end),col="red")+
    geom_point(data = filter(crisprs,seqid%in%cas_res$seqid),aes(x=start,y=end),col="blue")+
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10))+
    facet_wrap(.~seqid,scales = "free",ncol = 1)

  #尝试画一下，类似基因图
  library(gggenes)
  ?gggenes::geom_gene_arrow
  ggplot(example_genes, aes(xmin = start, xmax = end,
                            y = molecule, fill = gene)) +
    geom_gene_arrow() +
    facet_wrap(~ molecule, scales = "free")
}

