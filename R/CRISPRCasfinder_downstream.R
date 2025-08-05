#' Prepare the result files from CRISPRCasFinder
#'
#' @param input_folder the folder of CRISPRCasFinder result, contains GFF/,TSV/,result.json,...
#' @param output_folder output folder, default: NULL, don not output to file.
#' @param genome_name the genome_name
#' @param verbose verbose, default: TRUE
#'
#' @return crispr
#' @import dplyr
#' @export
#'
#' @examples
#' MAG_test <- system.file("extdata/MAG_test", package = "iCRISPR")
#' crispr <- pre_CCF_res(MAG_test)
#' @references https://crisprcas.i2bc.paris-saclay.fr/CrisprCasFinder/Index
pre_CCF_res <- function(input_folder, output_folder = NULL, genome_name = NULL, verbose = TRUE) {
  if (!dir.exists(input_folder)) stop("Can not find ", input_folder)
  if (!all(c("TSV", "GFF", "result.json") %in% list.files(input_folder))) stop("Maybe not a CRISPRCasFinder result.")

  if (is.null(genome_name)) genome_name <- basename(input_folder)

  # 有annotation前缀的gff文件是有crispr array的contigs
  crispr_gff_ls <- list.files(paste0(input_folder, "/GFF"), pattern = "annotation*") %>%
    gsub("annotation_", "", .)

  if (length(crispr_gff_ls) > 0) {
    crispr_res <- get_crispr(paste0(input_folder, "/TSV/Crisprs_REPORT.tsv"), genome_name = genome_name)
    check_crispr <- FALSE

    array_res <- res <- data.frame()
    for (i in crispr_gff_ls) {
      crispr_gff <- paste0(input_folder, "/GFF/", i)
      # tmp=get_spacer_fa(crispr_gff,genome_name=genome_name)
      # res=rbind(res,tmp)
      tmp <- get_array(crispr_gff, genome_name = genome_name)
      if (is.null(tmp)) check_crispr <- TRUE
      array_res <- rbind(array_res, tmp)
    }
    if (check_crispr) {
      crispr_res <- dplyr::filter(crispr_res, seqid %in% array_res$seqid)
    }
    if (nrow(array_res) < 2) {
      array_res <- spacer_res <- NULL
      if (verbose) pcutils::dabiao("No array found!")
    } else {
      spacer_res <- dplyr::filter(array_res, feature == "CRISPRspacer")
      tmp <- crispr_res %>%
        dplyr::mutate(long_id = paste0(genome, "@", CRISPR_id, "@", start, "-", end)) %>%
        dplyr::select(CRISPR_id, long_id)
      tmp1 <- dplyr::left_join(spacer_res, tmp, by = c("CRISPR_id" = "CRISPR_id")) %>%
        dplyr::group_by(CRISPR_id) %>%
        dplyr::mutate(sid = 1:n()) %>%
        dplyr::ungroup()
      spacer_res <- dplyr::mutate(tmp1, Spacer_id = paste0(long_id, "@spacer_", start, "_", abs(start - end) + 1, "_", sid))
      spacer_res <- dplyr::select(spacer_res, -long_id, -sid) %>% as.data.frame()
    }
  } else {
    array_res <- spacer_res <- crispr_res <- NULL
    if (verbose) pcutils::dabiao("No Crispr, Array found!")
  }

  rawCas_file <- paste0(input_folder, "/rawCas.fna")
  Cas_report_file <- paste0(input_folder, "/TSV/Cas_REPORT.tsv")
  if (file.exists(rawCas_file)) {
    cas_res <- get_cas(rawCas_file, Cas_report_file, genome_name = genome_name)
  } else {
    cas_res <- NULL
    if (verbose) pcutils::dabiao("No Cas found!")
  }

  # 导出包含cas的spacer
  # df=read_fasta(out_file1)
  # spacer_with_cas=df[grepl(paste0(cas_res$seqid,collapse = "|"),df$Sequence_ID),]
  # write_fasta(spacer_with_cas,out_file2)

  crispr <- structure(list(CRISPR = crispr_res, Cas = cas_res, Array = array_res, Spacer = spacer_res), class = "crispr")
  attributes(crispr)$basic_info <- list(
    genome_name = genome_name,
    n_crispr = length(unique(crispr$CRISPR$CRISPR_id)),
    n_cas = length(unique(crispr$Cas$Cas_id)),
    n_spacer = ifelse(is.null(nrow(crispr$Spacer)), 0, nrow(crispr$Spacer))
  )
  if (!is.null(output_folder)) {
    if (!output_folder == "") {
      output_folder <- paste0(output_folder, "/", genome_name)
      if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)
      out_file1 <- paste0(output_folder, "/", genome_name, "_spacer.fa")
      # 导出CRISPR信息
      if (!is.null(crispr_res)) utils::write.csv(crispr_res, file = paste0(output_folder, "/", genome_name, "_CRISPR_info.csv"), row.names = FALSE)
      # 导出Array信息
      if (!is.null(array_res)) utils::write.csv(array_res, file = paste0(output_folder, "/", genome_name, "_Array_info.csv"), row.names = FALSE)
      # 导出Spacer信息
      if (!is.null(spacer_res)) pcutils::write_fasta(dplyr::select(spacer_res, Spacer_id, sequence) %>% as.data.frame(), file_path = out_file1)
      if (!is.null(spacer_res)) utils::write.csv(spacer_res, file = paste0(output_folder, "/", genome_name, "_Spacer_info.csv"), row.names = FALSE)
      # 导出Cas信息
      if (!is.null(cas_res)) utils::write.csv(cas_res, file = paste0(output_folder, "/", genome_name, "_Cas_info.csv"), row.names = FALSE)
      if (verbose) pcutils::dabiao("All done! see ", output_folder, "/")
    }
  }
  if (verbose) pcutils::dabiao("All done!")
  return(crispr)
}

#' Prepare many genome result files from CRISPRCasFinder
#'
#' @param input_folder the folder of CRISPRCasFinder result, contains GFF/,TSV/,result.json,...
#' @param output_folder output folder, default: NULL, don not output to file.
#' @param threads threads, default: 1
#' @param verbose verbose
#'
#' @export
#' @return No value
multi_pre_CCF_res <- function(input_folder, output_folder = NULL, threads = 1, verbose = TRUE) {
  all_genome <- list.dirs(input_folder, recursive = FALSE)
  if (length(all_genome) == 0) {
    return(NULL)
  }
  all_genome <- all_genome[sapply(all_genome, \(i)all(c("TSV", "GFF", "result.json") %in% list.files(i)))]
  if (length(all_genome) == 0) {
    return(NULL)
  }
  pcutils::dabiao("Total ", length(all_genome), " genomes")
  pcutils::dabiao("Using ", threads, " threads")
  pcutils::dabiao("Start: ", date())
  # parallel
  reps <- length(all_genome)
  # main function
  loop <- function(i) {
    if (i %% 100 == 0) {
      pcutils::dabiao("Doing ", i, " genome")
    }

    tryCatch(
      {
        pre_CCF_res(
          input_folder = all_genome[i], output_folder = output_folder,
          verbose = FALSE
        )
      },
      error = function(e) {
        message("Error: ", all_genome[i])
        NULL
      }
    )
  }
  {
    if (threads > 1) {
      pcutils::lib_ps("foreach", "doSNOW", "snow", library = FALSE)
      if (verbose) {
        pb <- utils::txtProgressBar(max = reps, style = 3)
        opts <- list(progress = function(n) utils::setTxtProgressBar(pb, n))
      } else {
        opts <- NULL
      }
      cl <- snow::makeCluster(threads)
      doSNOW::registerDoSNOW(cl)
      res <- foreach::`%dopar%`(
        foreach::foreach(i = seq_len(reps), .options.snow = opts),
        loop(i)
      )
      snow::stopCluster(cl)
      gc()
    } else {
      res <- lapply(seq_len(reps), loop)
    }
  }
  # simplify method
  names(res) <- basename(all_genome)
  pcutils::dabiao("End: ", date())
  class(res) <- "multi_crispr"
  res
}

# deprecated，可以直接从array中生成所有spacer信息，不用从这里提取
# get_spacer_fa=function(crispr_gff,genome_name){
#   #用@分割，header要包括这些信息：spacer来源于哪条genome，哪条序列，是在该序列第几个crispr array 上，crispr array的起始和终止位点，再加一个spacer的起始位点，spacer的长度，以及是第几个spacer
#   res=data.frame(spacer_id=c(),spacer_seq=c())
#   gff1=utils::read.table(crispr_gff,col.names = c("seqid", "source", "feature", "start", "end", "score", "strand", "phase", "attributes"),comment.char = "#", header = FALSE, stringsAsFactors = FALSE,quote = "")
#   gff2=dplyr::filter(gff1,feature%in%c("CRISPR","CRISPRspacer"))%>%dplyr::select(1,3,4,5,9)
#   cid=0
#   for (r in 1:nrow(gff2)) {
#     row=gff2[r,]
#     if(row$feature=="CRISPR"){
#       cid=cid+1
#       sid=0
#       crispr1=row
#     }
#     else if(row$feature=="CRISPRspacer") {
#       sid=sid+1
#       attri=strsplit(row$attributes,";")[[1]]
#       spacer_id=paste0(genome_name,
#                        "@",row$seqid,
#                        "@CRISPR:",cid,
#                        "@",crispr1$start,"-",crispr1$end,
#                        "@",gsub("Name=","",attri[2]),"_",sid)
#       spacer_seq=gsub("sequence=","",attri[1])
#       # if(is.null(res))res=paste0(spacer_id,"\n",spacer_seq)
#       # else res=paste0(res,"\n",spacer_id,"\n",spacer_seq)
#       res=rbind(res,data.frame(spacer_id=spacer_id,spacer_seq=spacer_seq))
#     }
#   }
#   return(res)
# }


# 检查array有效性
check_array <- function(array) {
  crispr_pos1 <- dplyr::filter(array, feature == "LeftFLANK") %>% dplyr::pull(start)
  crispr_pos2 <- dplyr::filter(array, feature == "RightFLANK") %>% dplyr::pull(end)
  all(dplyr::between(range(c(array$start, array$end)), crispr_pos1, crispr_pos2))
}

get_array <- function(crispr_gff, genome_name) {
  gff <- pcutils::read.file(crispr_gff, all_yes = TRUE)
  # 部分结果有注释到crispr但是gff为空
  if (nrow(gff) < 2) {
    return(NULL)
  }

  # 必须做这个排序才能对应
  # gff=dplyr::arrange(gff,start)
  # 排序发现了bug，有的RightFLANK会跑到LeftFLANK前面，导致错位,因为有些spacer在array的外面？
  # 具体看output/dir_173/GCF_001390335.1_8939_3_24_genomic.fna/GFF/NZ_CPDV01000060.gff

  left <- which(gff$feature == "LeftFLANK")
  right <- which(gff$feature == "RightFLANK")
  if ((length(left) != length(right))) {
    warning("Something wrong with this gff file: ", crispr_gff)
    return(NULL)
  }
  # 新策略，先分开，再给id
  crispr_ls <- list()
  for (i in seq_along(left)) {
    crispr_ls[[i]] <- gff[left[i]:right[i], ]
  }

  if (!all(sapply(crispr_ls, check_array))) {
    warning("Something wrong with this gff file: ", crispr_gff, call. = FALSE)
    return(NULL)
  }

  ord <- rank(sapply(crispr_ls, \(i)i[1, "start"]))
  for (i in seq_along(left)) {
    crispr_ls[[i]]$CRISPR_id <- paste0("CRISPR:", ord[i])
  }
  gff <- do.call(rbind, crispr_ls)

  gff$CRISPR_id <- paste0(gff$seqid, "@", gff$CRISPR_id)
  gff <- dplyr::filter(gff, feature != "CRISPR")
  array_res <- gff[, c(1, 10, 3:5, 7)]
  array_res$sequence <- sub(".*sequence=([^;]+).*", "\\1", gff$attributes)
  array_res$strand <- ifelse(array_res$strand == "+", "Forward", ifelse(array_res$strand == "-", "Reverse", "Unknown"))
  return(data.frame(genome = genome_name, array_res))
}

# deprecated，不从这里提取
# get_cas=function(rawCas_file, Cas_report_file,genome_name){
#   rawCas_fasta=read_fasta(rawCas_file)%>%distinct(Sequence_ID,.keep_all = TRUE)
#   rawCas_fasta$Sequence_ID=lapply(rawCas_fasta$Sequence_ID,\(i)strsplit(i,split="\\|")[[1]][2])%>%do.call("c",.)
#   rawCas_fasta=rawCas_fasta%>%distinct(Sequence_ID,.keep_all = TRUE)
#
#   cas_report=readLines(Cas_report_file)
#   cas_report=cas_report[!(grepl("^[#-]",cas_report)|grepl("^$",cas_report)|grepl("\\( Unknown \\)",cas_report))]%>%
#     read.table(text = .,check.names = FALSE,
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

get_cas <- function(rawCas_file, Cas_report_file, genome_name) {
  all_cas_sys <- data.frame()
  cas_report <- readLines(Cas_report_file)
  cas_sys <- c()
  flag <- FALSE
  flag2 <- FALSE
  contig <- ""

  for (i in cas_report) {
    # System开始部分
    if (grepl("### System:", i)) {
      sys_type <- sub("### System: (.*) \\(.*\\)", "\\1", i)
      flag <- TRUE
    }
    # System结束部分
    else if (grepl("####Summary system", i)) {
      flag <- FALSE
      # 有些Sequence_ID内带#号，恶心坏了
      cas_sys <- cas_sys[!grepl("^#", cas_sys)]
      a <- utils::read.table(
        text = cas_sys, check.names = FALSE, comment.char = "", quote = "", sep = "\t",
        col.names = c("Sequence_ID", "Cas-type/subtype", "Gene status", "System", "Type", "Begin", "End", "Strand", "Other_information")
      )
      a <- a %>%
        dplyr::select("Sequence_ID", "Type", "System", "Cas-type/subtype", "Begin", "End", "Strand") %>%
        dplyr::distinct()

      contig1 <- sub(".*sequenceID=(\\S+)}.*", "\\1", i)
      if (contig1 != contig) {
        cas_id <- 1
      } else {
        cas_id <- cas_id + 1
      }
      contig <- contig1

      a$Cas_id <- paste0("CAS:", cas_id)
      a$Type <- sys_type
      cas_sys <- c()
      all_cas_sys <- rbind(all_cas_sys, a)
    }
    # System开始时读入新行
    if (flag) cas_sys <- c(cas_sys, i)
  }

  rawCas_fasta <- pcutils::read_fasta(rawCas_file) %>% dplyr::distinct(Sequence_ID, .keep_all = TRUE)
  rawCas_fasta$Sequence_ID <- lapply(rawCas_fasta$Sequence_ID, \(i)strsplit(i, split = "\\|")[[1]][2]) %>% do.call("c", .)
  rawCas_fasta <- rawCas_fasta %>% dplyr::distinct(Sequence_ID, .keep_all = TRUE)
  cas_res <- dplyr::left_join(all_cas_sys, rawCas_fasta, by = c("Sequence_ID" = "Sequence_ID"))
  cas_res$feature <- "Cas"
  colnames(cas_res) <- c("seqid", "type", "subtype", "protein", "start", "end", "strand", "Cas_id", "sequence", "feature")
  cas_res$strand <- ifelse(cas_res$strand == "+", "Forward", ifelse(cas_res$strand == "-", "Reverse", "Unknown"))
  cas_res$seqid <- gsub("_\\d+$", "", cas_res$seqid)

  cas_res$subtype <- gsub("CAS-", "", cas_res$type)
  cas_res$type <- ifelse(grepl("Type", cas_res$subtype), substr(cas_res$subtype, 1, nchar(cas_res$subtype) - 1), cas_res$subtype)

  cas_res$Cas_id <- paste0(cas_res$seqid, "@", cas_res$Cas_id)
  cas_res <- cbind(genome = genome_name, cas_res[, c("seqid", "Cas_id", "feature", "type", "subtype", "protein", "start", "end", "strand", "sequence")])
  cas_res
}

get_crispr <- function(Crisprs_REPORT, genome_name = genome_name) {
  # Crisprs_REPORT="pre_CCF_res_out/GCA_005025685.1_PDT000277779.2_genomic.fna/TSV/Crisprs_REPORT.tsv"
  crisprs <- utils::read.table(Crisprs_REPORT, sep = "\t", header = TRUE, check.names = FALSE, comment.char = "", stringsAsFactors = FALSE, quote = "")
  crisprs <- crisprs[, c(
    "Sequence", "CRISPR_Id", "Strain", "CRISPR_Start", "CRISPR_End",
    "Potential_Orientation (AT%)", "Consensus_Repeat", "Spacers_Nb", "Evidence_Level"
  )]
  crisprs$Strain <- "CRISPR"
  crisprs$CRISPR_Id <- paste0(crisprs$Sequence, "@CRISPR:", sub(".*_(\\d+)$", "\\1", crisprs$CRISPR_Id))
  colnames(crisprs) <- c("seqid", "CRISPR_id", "feature", "start", "end", "strand", "consensus_repeat", "spacer_number", "evidence_level")
  data.frame(genome = genome_name, crisprs)
}

#' Plot a crispr-cas system
#'
#' @param crispr crispr result from `pre_CCF_res()`
#' @param genome the genome name
#' @param contig the contig name
#' @param array plot the array?
#' @param cas plot the cas?
#'
#' @import ggplot2
#' @return ggplot
#' @export
#'
#' @examples
#' data(crispr)
#' (p <- plot_crispr(crispr,
#'   genome = "MAG_test",
#'   contig = "AAB-S01R1_k55_9399631_flag=0_multi=9.8751_len=26518"
#' ))
#' p + xlim_crispr(24000, 27000)
plot_crispr <- function(crispr, contig = NULL, genome = NULL, array = TRUE, cas = TRUE) {
  if (!inherits(crispr, "crispr")) {
    return(NULL)
  }
  with_cas <- !is.null(crispr$Cas)
  with_array <- !is.null(crispr$CRISPR)
  if (!(with_cas | with_array)) {
    return(NULL)
  }
  if (is.null(contig)) {
    message("Please set the contig as one seqid in crispr!")
    if (with_array) contig <- crispr$CRISPR$seqid[1]
    if (with_cas) contig <- crispr$Cas$seqid[1]
    message("Use ", contig, " as contig")
  }

  # 尝试画一下，类似基因图
  pcutils::lib_ps("ggnewscale", "gggenes", library = FALSE)
  cas_res <- data.frame()
  array_res <- data.frame()
  if (with_cas & cas) {
    cas_res <- crispr$Cas %>% dplyr::filter(seqid == contig)
    cas_res <- dplyr::mutate(cas_res, protein = gsub("_.*", "", protein))
  }
  # crispr_res=crispr$CRISPR%>%dplyr::filter(seqid==contig)
  if (with_array & array) {
    array_res <- crispr$Array %>% dplyr::filter(seqid == contig)
    array_res$feature <- factor(array_res$feature, levels = c("LeftFLANK", "CRISPRdr", "CRISPRspacer", "RightFLANK")) %>% droplevels()
  }
  if (is.null(genome)) genome <- attributes(crispr)$basic_info$genome_name
  sub_title <- ""

  p <- ggplot2::ggplot()
  if (nrow(array_res) > 0) {
    spacer_n <- array_res %>%
      dplyr::filter(feature == "CRISPRspacer") %>%
      nrow()
    sub_title <- paste0(sub_title, "Spacer number: ", spacer_n)

    p <- p +
      # crispr系统
      gggenes::geom_gene_arrow(
        data = array_res, aes(xmin = start, xmax = end, y = seqid, fill = feature, forward = (strand != "Reverse")),
        color = NA, arrowhead_width = grid::unit(0, "mm"), arrowhead_height = grid::unit(0, "mm")
      ) +
      scale_fill_manual(values = setNames(
        c("#5950FA", "#190861", "#78c679", "#fb8072"),
        c("LeftFLANK", "CRISPRdr", "CRISPRspacer", "RightFLANK")
      )) +
      ggnewscale::new_scale_fill()
  }
  if (nrow(cas_res) > 0) {
    cas_n <- dplyr::distinct(cas_res, Cas_id, subtype)
    sub_title <- paste0(sub_title, "Cas-system number: ", nrow(cas_n), "; Subtype: ", paste0(cas_n$subtype, collapse = "/"), "; ")

    p <- p +
      # cas基因
      gggenes::geom_gene_arrow(
        data = cas_res, aes(xmin = start, xmax = end, y = seqid, fill = protein, forward = (strand != "Reverse")),
        arrowhead_width = grid::unit(5, "mm"), arrowhead_height = grid::unit(5, "mm"), arrow_body_height = grid::unit(4, "mm")
      ) +
      gggenes::geom_gene_label(data = cas_res, aes(xmin = start, xmax = end, y = seqid, fill = protein, label = protein), min.size = 0) +
      scale_fill_manual(values = pcutils::get_cols(length(unique(cas_res$protein))))
  }
  if (!((nrow(cas_res) > 0) | (nrow(array_res) > 0))) {
    message("No cas or crispr in this contig: ", contig)
  } else {
    p <- p + coord_fixed(diff(range(c(cas_res$start, array_res$start, cas_res$end, array_res$end))) / 8)
  }
  p <- p + gggenes::theme_genes() +
    labs(y = NULL, title = paste0("Genome: ", genome, "; Contig: ", contig), subtitle = sub_title) +
    theme(
      legend.position = "top",
      axis.line.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  p
}

# plot_crispr=function(crispr,CRISPR_id=NULL,genome=NULL,array=TRUE,cas=TRUE){
#   if(!inherits(crispr,"crispr"))return(NULL)
#   if(is.null(crispr$CRISPR))return(NULL)
#   if(is.null(CRISPR_id)){
#     message("Please set the CRISPR_id as one CRISPR_id in crispr!")
#     contig=crispr$CRISPR$CRISPR_id[1]
#     message("Use ",contig," as CRISPR_id")
#   }
#   #尝试画一下，类似基因图
#   array_res=data.frame()
#   array_res=crispr$Array%>%dplyr::filter(seqid==contig)
#   array_res$feature=factor(array_res$feature,levels = c("LeftFLANK","CRISPRdr","CRISPRspacer","RightFLANK"))%>%droplevels()
#   if(is.null(genome))genome=attributes(crispr)$basic_info$genome_name
#   sub_title=""
#
# }

#' xlim for plot_crispr
#'
#' @param ... numeric
#'
#' @return ggproto object
#' @export
#' @rdname plot_crispr
xlim_crispr <- function(...) {
  ggplot2::coord_cartesian(xlim = c(...))
}

#' Get spacer fasta
#'
#' @param crispr crispr
#' @param evidence_level filter the evidence_level (1~4), if NULL, then no filter
#' @param cas filter with cas system or not (TRUE or FALSE), if NULL, then no filter
#'
#' @export
#' @return fasta
#' @examples
#' data(crispr)
#' get_spacer_fa(crispr, evidence_level = 4, cas = TRUE)
get_spacer_fa <- function(crispr, evidence_level = NULL, cas = NULL) {
  if (inherits(crispr, "multi_crispr")) {
    aaa <- lapply(crispr, get_spacer_fa, evidence_level = evidence_level, cas = cas)
    bbb <- do.call(rbind, aaa)
    rownames(bbb) <- NULL
    return(bbb)
  }
  if (!inherits(crispr, "crispr")) {
    return(NULL)
  }
  if (is.null(crispr$Spacer)) {
    return(NULL)
  }

  if (!is.null(cas)) {
    if (cas & is.null(crispr$Cas)) {
      return(NULL)
    }
    if (!cas & !is.null(crispr$Cas)) {
      return(NULL)
    }
  }

  if (!is.null(evidence_level)) {
    # 旧版本生成名称为evidence_Level，已改为evidence_Level
    if ("evidence_Level" %in% colnames(crispr$CRISPR)) crispr$CRISPR <- dplyr::rename(crispr$CRISPR, evidence_level = "evidence_Level")
    filter_crispr <- dplyr::filter(crispr$CRISPR, evidence_level %in% !!evidence_level) %>% dplyr::pull(CRISPR_id)
    filter_spacer <- crispr$Spacer %>%
      dplyr::filter(CRISPR_id %in% filter_crispr) %>%
      dplyr::select(Spacer_id, sequence)
  } else {
    filter_spacer <- crispr$Spacer %>% dplyr::select(Spacer_id, sequence)
  }
  if (nrow(filter_spacer) == 0) {
    return(NULL)
  }
  return(filter_spacer)
}

#' Get consensus repeat fasta
#'
#' @param crispr crispr
#' @param evidence_level filter the evidence_level (1~4), if NULL, then no filter
#' @param cas filter with cas system or not (TRUE or FALSE), if NULL, then no filter
#'
#' @export
#' @return fasta
#' @examples
#' data(crispr)
#' get_consensus_rep_fa(crispr)
get_consensus_rep_fa <- function(crispr, evidence_level = NULL, cas = NULL) {
  if (inherits(crispr, "multi_crispr")) {
    aaa <- lapply(crispr, get_consensus_rep_fa, evidence_level = evidence_level, cas = cas)
    bbb <- do.call(rbind, aaa)
    rownames(bbb) <- NULL
    return(bbb)
  }
  if (!inherits(crispr, "crispr")) {
    return(NULL)
  }
  if (is.null(crispr$CRISPR)) {
    return(NULL)
  }

  if (!is.null(cas)) {
    if (cas & is.null(crispr$Cas)) {
      return(NULL)
    }
    if (!cas & !is.null(crispr$Cas)) {
      return(NULL)
    }
  }

  if (!is.null(evidence_level)) {
    # 旧版本生成名称为evidence_Level，已改为evidence_Level
    if ("evidence_Level" %in% colnames(crispr$CRISPR)) crispr$CRISPR <- dplyr::rename(crispr$CRISPR, evidence_level = "evidence_Level")
    filter_consensus_rep <- dplyr::filter(crispr$CRISPR, evidence_level == !!evidence_level) %>% dplyr::select(CRISPR_id, consensus_repeat)
  } else {
    filter_consensus_rep <- crispr$CRISPR %>% dplyr::select(CRISPR_id, consensus_repeat)
  }
  if (nrow(filter_consensus_rep) == 0) {
    return(NULL)
  }
  return(filter_consensus_rep)
}

update_crispr <- function(crispr) {
  # 有一些版本更替的变化放在这里
  # 比如Cas的TypeU，应该让type也变成TypeU
  if (inherits(crispr, "multi_crispr")) {
    aaa <- lapply(crispr, update_crispr)
    class(aaa) <- c("multi_crispr")
    return(aaa)
  }
  if (!inherits(crispr, "crispr")) {
    return(NULL)
  }

  if (is.null(crispr$Cas)) {
    return(crispr)
  }
  check_u <- crispr$Cas$subtype == "TypeU"
  if (any(check_u)) {
    crispr$Cas$type[which(check_u)] <- "TypeU"
  }
  return(crispr)
}

#' Combine some kinds of objects in `iCRISPR`
#'
#' @param ... additional
#'
#' @return same object
#' @export
#'
combine <- function(...) {
  UseMethod("combine")
}

#' Combine some crispr or multi_crispr object
#'
#' @param ... crispr or multi_crispr object
#'
#' @return multi_crispr object
#' @exportS3Method
#' @method combine crispr
#'
combine.crispr <- function(...) {
  combine_crispr(...)
}

#' Combine some crispr or multi_crispr object
#'
#' @param ... crispr or multi_crispr object
#'
#' @return multi_crispr object
#' @exportS3Method
#' @method combine multi_crispr
#'
combine.multi_crispr <- function(...) {
  combine_crispr(...)
}

combine_crispr <- function(...) {
  res <- list()
  all_c <- list(...)
  for (crispr in all_c) {
    if (inherits(crispr, "multi_crispr")) {
      res <- append(res, crispr)
    }
    if (inherits(crispr, "crispr")) {
      res[[attributes(crispr)$basic_info$genome_name]] <- crispr
    }
  }
  class(res) <- "multi_crispr"
  return(res)
}

#' Select part of multi_crispr.
#'
#' @param .data multi_crispr
#' @param ... additional
#'
#' @return multi_crispr
#'
#' @export
#'
#' @examples
#' data(multi_crispr)
#' select_crispr(multi_crispr, 1:10)
select_crispr <- function(.data, ...) {
  if (!inherits(.data, "multi_crispr")) warning("not multi_crispr!")
  selected <- .data[...]
  class(selected) <- class(.data)
  return(selected)
}
