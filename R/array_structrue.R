#' Get Arrays from crispr.
#'
#' @param crispr crispr object
#'
#' @return Arrays
#' @export
#'
#' @examples
#' data(crispr)
#' Arrays=get_Arrays(crispr)
get_Arrays=function(crispr){
  crispr_id=unique(crispr$Array$CRISPR_id)

  Arrays=lapply(crispr_id,\(i){
    tmp=crispr$CRISPR%>%dplyr::filter(CRISPR_id==i)
    res=as.list(tmp)
    spacers=crispr$Spacer%>%dplyr::filter(CRISPR_id==i)%>%dplyr::pull(sequence)
    res=append(res,list(spacers=spacers))
    class(res)="Array"
    res
  })
  names(Arrays)=crispr_id
  Arrays
}

#' Print Array object
#'
#' @method print Array
#' @param ... add
#' @param x Array object
#'
#' @exportS3Method
print.Array=function(x,...){
  pcutils::dabiao("CRISPR_id: ",x$CRISPR_id,print = TRUE)
  cat("With",x$spacer_number,"spacers; Evidence level=",x$evidence_level,"\n")
  cat("With consensus_repeat:",x$consensus_repeat,"\n")
}

get_whole_seq=function(x_pattern){
  aligned_pattern <- Biostrings::aligned(x_pattern)[[1L]]
  ans <- aligned_pattern
  original_pattern <- x_pattern@unaligned[[1L]]
  start1 <- Biostrings::start(x_pattern@range)
  if (start1 > 1L) {
    prefix1 <- Biostrings::subseq(original_pattern, end = start1 - 1L)
    ans <- c(prefix1, ans)
  }
  end1 <- Biostrings::end(x_pattern@range)
  if (end1 < length(original_pattern)) {
    suffix1 <- Biostrings::subseq(original_pattern, start = end1 +1L)
    ans <- c(ans, suffix1)
  }
  ans
}

#' Calculate Sequence Identity of two sequence
#'
#' @description
#' This function calculates the sequence identity between two sequences.
#'
#' @param s1 First sequence
#' @param s2 Second sequence
#' @param xmode consider the reverse string (all), complement and reverse complement string (DNA)
#'
#' @return Sequence identity as a decimal value between 0 and 1
#' @export
#' @examples
#' calculate_identity("ATCGTACG","ATCGTAGC")
#' calculate_identity("ATCGTACG","ATCGTAGC",xmode=TRUE)
calculate_identity <- function(s1, s2,xmode=FALSE) {
  pcutils::lib_ps("Biostrings",library = FALSE)
  if(xmode){
    s1s=c(s1,Biostrings::reverse(s1))
    flag=tryCatch(expr = {dna_s1=Biostrings::DNAString(s1);TRUE},error=function(e){FALSE})
    if(flag){
      s1s=c(s1s,Biostrings::complement(Biostrings::DNAString(s1)))
      s1s=c(s1s,Biostrings::reverseComplement(Biostrings::DNAString(s1)))
    }
    identity=sapply(s1s, \(i)calculate_identity(i,s2))
    return(identity)
  }
  s1=as.character(s1)
  s2=as.character(s2)
  # Get all unique characters from both strings
  all_chars <- unique(c(strsplit(s1, "")[[1]], strsplit(s2, "")[[1]]))
  # Create a substitution matrix with diagonal elements as 1
  mat <- matrix(diag(length(all_chars)), nrow = length(all_chars), ncol = length(all_chars), dimnames = list(all_chars, all_chars))

  # Calculate identity
  globalAlign <- Biostrings::pairwiseAlignment(s1, s2, substitutionMatrix = mat, gapOpening = 0, gapExtension = 0)
  #提供的函数pid只加了内部的gap，不合理
  #identity=Biostrings::pid(globalAlign)

  #ans=Biostrings:::get_aligned_pattern(globalAlign@pattern,globalAlign@subject)

  ans=get_whole_seq(globalAlign@pattern)
  whole_length=Biostrings::nchar(ans)
  identity=Biostrings::nmatch(globalAlign)/whole_length
  return(identity)
}

#' Compare two crispr arrays (whatever arrays indeed)
#'
#' @param array1 a vector contains spacers (with order) as reference.
#' @param array2 a vector contains spacers (with order)
#'
#' @return array_comparison object
#' @export
#'
#' @examples
#' array_test=random_seq(5)[,2]
#' res=compare_array(array1=array_test,array2=array_test[c(5,1:3)])
#' plot.array_comparison(res)
#' align_array(res)->align_res
#' plot.array_comparison(align_res)
compare_array=function(array1,array2){
  array1_id=names(array1)
  array2_id=names(array2)
  if(!(is.null(array1_id)|is.null(array2_id))){
    if(length(intersect(array1_id,array2_id))>0)stop("check name of array1 and array2")
  } else {
    array1_id=paste0("A",seq_along(array1))
    array2_id=paste0("B",seq_along(array2))
  }

  array_spacer=rbind(data.frame(id=array1_id,spacer=array1,array="array1",position=seq_along(array1)),
                     data.frame(id=array2_id,spacer=array2,array="array2",position=seq_along(array2)))

  compare=expand.grid(array1_id,array2_id)
  colnames(compare)=paste0("id",1:2)
  compare=left_join(compare,array_spacer,by=c(id1="id"))
  compare=left_join(compare,array_spacer,by=c(id2="id"),suffix = c("1", "2"))

  #compare$identity=apply(compare,1,\(i)calculate_identity(i[3],i[6]))
  compare$identity=apply(compare,1,\(i)max(calculate_identity(i[3],i[6],xmode = TRUE)))

  compare=compare%>%mutate(Identity_level = cut(identity, breaks = c(-Inf, 0.8, 0.9, Inf),
                                    labels = c("< 0.8", "0.8 - 0.9", ">= 0.9")))

  #把有link的spacer归为一组
  links=filter(compare,identity>0.9)
  if(nrow(links)>0){
    pcutils::lib_ps("igraph",library = FALSE)
    spacer_net=igraph::graph_from_data_frame(links,vertices = array_spacer)
    #把有link的spacer归为一组
    spacer_membership=igraph::clusters(spacer_net)%>%igraph::membership()
    array_spacer$membership=paste0("M",spacer_membership[array_spacer$id])
  }
  else array_spacer$membership=paste0("M",seq_len(nrow(array_spacer)))

  res=list(array_spacer=array_spacer,compare=compare)
  class(res)="array_comparison"
  return(res)
}


#' Plot array_comparison
#'
#' @param x array_comparison
#' @param ... add
#' @param linewidth linewidth
#'
#' @return ggplot
#' @exportS3Method
#' @method plot array_comparison
#'
#' @rdname compare_array
plot.array_comparison=function(x,linewidth=2,...){
  compare=x$compare
  array_spacer=x$array_spacer
  compare$Identity_level=factor(compare$Identity_level,levels = rev(c("< 0.8", "0.8 - 0.9", ">= 0.9")))
  links=filter(compare,identity>0.9)
  p=ggplot()
  if(nrow(links)>0){
    p=p+geom_segment(data = dplyr::filter(compare,identity>=0.8),
                            mapping = aes(x=position1,xend=position2,y=array1,yend=array2,col=Identity_level),
                     linewidth=linewidth,alpha=0.8)+
      scale_color_manual(values = setNames(c("#1f78b4", "#ffed6f","#fb8072"),c("< 0.8", "0.8 - 0.9", ">= 0.9")))
  }

  if(any(duplicated(array_spacer$membership))){
    p=p+
      geom_tile(data = array_spacer,mapping = aes(x=position,y=array,fill=membership),width=0.55,height=0.25)+
      scale_fill_manual(values = pcutils::get_cols(length(unique(array_spacer$membership))),guide="none")
    }
  else p=p+geom_tile(data = array_spacer,mapping = aes(x=position,y=array),width=0.55,height=0.25,fill="#78c679")

  if(!is.null(x$array_align_res)){
    sv_res=summary_SV(x$array_align_res)
    title=paste0("Length: ",sv_res$Length,"; Match: ",sv_res$Match,"; Translocation: ",
                 sv_res$Translocation,"\nDeletion: ",sv_res$Deletion,"; Insertion: ",sv_res$Insertion)
    p=p+labs(title = title)
  }
  if(max(array_spacer$position)>20)p=p+coord_fixed(ratio = 1.5)
  else if(max(array_spacer$position)>40)p=p+coord_fixed(ratio = 2)
  else p=p+coord_fixed(ratio = 1)
  p+theme_bw()+labs(x="Position",y=NULL)
}

#' Align array_comparison
#'
#' @param array_comparison array_comparison object
#'
#' @return aligned array_comparison object
#' @export
#'
#' @rdname compare_array
align_array=function(array_comparison){
  compare=array_comparison$compare
  array_spacer=array_comparison$array_spacer

  #把有link的spacer归为一组，抽象成字符进行序列比对
  links=filter(array_comparison$compare,identity>0.9)
  if(nrow(links)>0){
    #抽象成字符进行序列比对
    seq1=array_comparison$array_spacer%>%filter(array=="array1")
    seq1=setNames(seq1$membership,seq1$id)
    seq2=array_comparison$array_spacer%>%filter(array=="array2")
    seq2=setNames(seq2$membership,seq2$id)

    array_align_res=Vector_alignment(seq1,seq2)
    array1=array_align_res$align1
    array2=array_align_res$align2
    array1_id=names(array1)
    array2_id=names(array2)

    align_array_spacer=rbind(data.frame(id=array1_id,spacer=array1,array="array1",position=seq_along(array1)),
                       data.frame(id=array2_id,spacer=array2,array="array2",position=seq_along(array2)))

    rownames(align_array_spacer)=NULL
    new_position=setNames(align_array_spacer$position,align_array_spacer$id)
    array_spacer$position=new_position[array_spacer$id]
    compare$position1=new_position[compare$id1]
    compare$position2=new_position[compare$id2]

    res=list(array_spacer=array_spacer,compare=compare,array_align_res=array_align_res)
    attributes(res)$align=TRUE
    class(res)="array_comparison"
    return(res)
  }
  else {
    message("no link!")
    return(res)
  }
}

#' Dynamic programming, a little slow but can used for multi-characteristics.
#'
#' @param seq1 vector
#' @param seq2 vector
#' @param mat substitution matrix
#' @param gap gap penalty
#' @param print_score_mat print score matrix
#'
#' @return alignment result
#' @export
#'
#' @examples
#' align_res=Vector_alignment(c("1"="as","2"="bb","3"="cc"),c("bb","cc","as","as"))
#' align_res
#' summary_SV(align_res)
#' Vector_alignment("ATCGTACG","ATCGTAGC",print_score_mat=TRUE)
Vector_alignment <- function(seq1, seq2, mat=NULL, gap =0,print_score_mat=FALSE) {
  #单元素就当作字符比较
  if(length(seq1)*length(seq2)==1){
    seq1=strsplit(seq1,"")[[1]]
    seq2=strsplit(seq2,"")[[1]]
  }
  # Get all unique characters from both strings
  all_chars <- unique(c(seq1,seq2))
  if(is.null(mat)){
    # Create a substitution matrix with diagonal elements as 1
    mat <- matrix(diag(length(all_chars)), nrow = length(all_chars), ncol = length(all_chars), dimnames = list(all_chars, all_chars))
  }
  else {
    if(!(all(all_chars%in%colnames(mat))&all(all_chars%in%rownames(mat))))stop("Check your substitution matrix")
  }

  n <- length(seq1)
  m <- length(seq2)

  # Create a matrix to store the scores
  score <- matrix(0, nrow = n + 1, ncol = m + 1)
  rownames(score)=c(0,seq1)
  colnames(score)=c(0,seq2)
  # Initialize the first row and first column with gap penalties
  score[1, ] <- seq(from = 0, to = -(m * abs(gap)), by = (gap))
  score[, 1] <- seq(from = 0, to = -(n * abs(gap)), by = (gap))

  # Fill in the matrix
  for (i in 2:(n + 1)) {
    for (j in 2:(m + 1)) {
      # Calculate the score for a match/mismatch
      # match_score <- score[i - 1, j - 1] + (ifelse(seq1[i-1] == seq2[j-1], match, mismatch))

      match_score <- score[i - 1, j - 1] + mat[seq1[i-1],seq2[j-1]]
      # Calculate the score for a gap in sequence 1
      gap1_score <- score[i - 1, j] + gap

      # Calculate the score for a gap in sequence 2
      gap2_score <- score[i, j - 1] + gap

      # Take the maximum score
      score[i, j] <- max(match_score, gap1_score, gap2_score)
    }
  }

  # Trace back to find the alignment
  align1 <- c()
  align2 <- c()
  i <- n + 1
  j <- m + 1
  if(print_score_mat){
    cat(strrep("=",100),"\n")
    print(score)
    cat(strrep("=",100),"\n")
  }
  f_score=score[i, j]

  while (i > 1 && j > 1) {
    current_score <- score[i, j]
    diag_score <- score[i - 1, j - 1]
    left_score <- score[i, j - 1]
    up_score <- score[i - 1, j]

    if (current_score == diag_score + mat[seq1[i-1],seq2[j-1]]) {
      align1 <- c(seq1[i-1], align1)
      align2 <- c(seq2[j-1], align2)
      i <- i - 1
      j <- j - 1
    } else if (current_score == up_score + gap) {
      align1 <- c(seq1[i-1], align1)
      align2 <- c("-", align2)
      i <- i - 1
    } else {
      align1 <- c("-", align1)
      align2 <- c(seq2[j-1], align2)
      j <- j - 1
    }
  }

  # If there are remaining characters in sequence 1, add them as gaps
  while (i > 1) {
    align1 <- c(seq1[i-1], align1)
    align2 <- c("-", align2)
    i <- i - 1
  }

  # If there are remaining characters in sequence 2, add them as gaps
  while (j > 1) {
    align1 <- c("-", align1)
    align2 <- c(seq2[j-1], align2)
    j <- j - 1
  }

  return(list(align1 = align1, align2 = align2,score=f_score))
}


#' Summary SV of an alignment
#'
#' @param align_res result from Vector_alignment
#'
#' @return list
#' @export
#'
#' @rdname Vector_alignment
summary_SV=function(align_res){
  align1=align_res$align1
  align2=align_res$align2
  seq1=align1[align1!="-"]
  seq2=align2[align2!="-"]
  Match=Deletion=Translocation=Insertion=0
  #count match
  Mp=(align1==align2)
  Match=sum(Mp)
  #resume
  align1=align1[!Mp]
  align1=align1[align1!="-"]
  align2=align2[!Mp]
  align2=align2[align2!="-"]
  for (i in seq_along(align2)) {
    if(align2[i]%in%align1){
      Translocation=Translocation+1
      align1=align1[-which(align1==align2[i])[1]]
    }
    else {
      Insertion=Insertion+1
    }
  }
  Deletion=length(align1)

  #定义一下两者结构变异的程度
  #match加一分，translocation0.5，insertion和deletion0分
  #总分为alignment后的长度。
  Length=length(align_res$align1)
  SV_index=(sum(Match)+0.5*sum(Translocation))/Length

  return(list(align1=align_res$align1,align2=align_res$align2,SV_index=SV_index,Length=Length,
              Match=Match,Deletion=Deletion,Translocation=Translocation,Insertion=Insertion))
}


#' Install MAFFT
#'
#' This function downloads and installs MAFFT at the specified path.
#'
#' @param mafft.zip your download mafft.zip file from https://mafft.cbrc.jp/alignment/software/
#' @param install_dir The directory to install MAFFT.
#'
#' @return None
#' @export
install_mafft <- function(mafft.zip=NULL,install_dir=tools::R_user_dir("iCRISPR")) {
  if (os <- tolower(Sys.info()[["sysname"]])!= "darwin")stop("Only for macos, please try to install mafft yourself from https://mafft.cbrc.jp/alignment/software/")

  install_dir=normalizePath(install_dir)
  # Create the installation directory if it does not exist
  dir.create(install_dir, recursive = TRUE, showWarnings = FALSE)

  # Set the full path to the installation directory
  install_path <- file.path(install_dir, "mafft.zip")

  url="https://mafft.cbrc.jp/alignment/software/mafft-7.490-mac.zip"
  if(is.null(mafft.zip)){
    # Download the file
    ori_time=getOption("timeout")
    on.exit(options(timeout = ori_time))

    options(timeout = 60)
    tryCatch(expr = {
      # Download MAFFT zip file
      download.file(url, destfile = install_path)
    },error=function(e){
      stop("Try download yourself from https://mafft.cbrc.jp/alignment/software/")
    })
  }
  else {
    if(file.exists(mafft.zip)&grepl("mafft.*\\.zip",mafft.zip)){
      filename <- basename(mafft.zip)
      install_path <- file.path(install_dir, filename)
      file.copy(mafft.zip,install_path)
    }
    else stop("Wrong file: ",mafft.zip)
  }

  # Unzip the downloaded file
  unzip(install_path)

  # Move the extracted folder to the installation path
  file.rename("mafft-mac", file.path(install_dir,"mafft-mac"))

  # Clean up temporary files
  file.remove(install_path)

  set_config("MAFFT",paste0(install_dir,"/mafft-mac/mafft.bat"))
  cat("MAFFT installed successfully at", install_path, "\n")
}
