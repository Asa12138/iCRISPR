get_crispr_info=function(crispr){
  if(is.null(crispr)) return(NULL)
  if (inherits(crispr, "multi_crispr")) {
    aaa <- lapply(crispr, get_crispr_info)
    bbb <- do.call(rbind, aaa)
    rownames(bbb) <- NULL
    return(bbb)
  }
  data.frame(attributes(crispr)$basic_info,
             n_cas_protein=ifelse(is.null(crispr$Cas),0,nrow(crispr$Cas)),
             n_crispr_level4=sum(crispr$CRISPR$evidence_level==4)
  )
}


#' Summary Cas-system type
#'
#' @param crispr crispr object
#' @param each_genome summary each genome information in multi_crispr?
#' @export
#' @return cas_type object
#' @examples
#' data(multi_crispr)
#' cas_type_res <- summary_cas_type(multi_crispr)
#' plot(cas_type_res)
summary_cas_type <- function(crispr, each_genome = FALSE) {
  if (inherits(crispr, "multi_crispr")) {
    aaa <- lapply(crispr, summary_cas_type)
    bbb <- do.call(rbind, aaa)
    if (is.null(bbb)) {
      return(NULL)
    }
    rownames(bbb) <- NULL

    if (!each_genome) {
      bbb <- dplyr::group_by(bbb[, c("type", "subtype", "n")], type, subtype) %>%
        dplyr::summarise(genome_num = dplyr::n(), n = sum(n)) %>%
        as.data.frame()
    }

    class(bbb) <- c("cas_type", class(bbb))
    return(bbb)
  }
  if (!inherits(crispr, "crispr")) {
    return(NULL)
  }
  cas_info <- crispr$Cas
  if (is.null(cas_info)) {
    return(NULL)
  }
  cas_type <- cas_info[, c("genome", "Cas_id", "type", "subtype")] %>%
    dplyr::distinct() %>%
    dplyr::count(genome, type, subtype) %>%
    as.data.frame()
  class(cas_type) <- c("cas_type", class(cas_type))
  return(cas_type)
}

#' Combine some cas_type object
#'
#' @param ... cas_type object
#' @param each_genome each_genome?
#'
#' @return cas_type object
#' @exportS3Method
#' @method combine cas_type
#'
combine.cas_type <- function(..., each_genome = TRUE) {
  all_c <- list(...)
  if (!all(sapply(all_c, inherits, what = "cas_type"))) stop("all input should be cas_type object")
  res <- lapply(all_c, \(i){
    if ("genome" %in% colnames(i)) {
      bbb <- dplyr::group_by(i[, c("type", "subtype", "n")], type, subtype) %>%
        dplyr::summarise(genome_num = dplyr::n(), n = sum(n)) %>%
        as.data.frame()
      return(bbb)
    } else {
      return(i)
    }
  })
  all_res <- do.call(rbind, res)
  ccc <- dplyr::group_by(all_res, type, subtype) %>%
    dplyr::summarise(genome_num = sum(genome_num), n = sum(n)) %>%
    as.data.frame()
  class(ccc) <- c("cas_type", "data.frame")
  ccc
}


#' Plot the Cas system type and subtype
#' @method plot cas_type
#'
#' @param x cas_type object
#' @param mode 1~3
#' @param ... additional cas_type object
#' @param plot_param plot parameters
#'
#' @return ggplot
#' @exportS3Method
#' @rdname summary_cas_type
plot.cas_type <- function(x, ..., mode = 1, plot_param = list()) {
  cas_type2 <- list(...)
  if (length(cas_type2) > 0) {
    type <- ifelse(mode == 1, "subtype", "type")
    cas_type <- lapply(
      append(list(x), cas_type2),
      \(x)(dplyr::group_by(x[, c(type, "n")], get(type)) %>%
        dplyr::summarise(n = sum(n)) %>% as.data.frame())
    )

    cas_type <- lapply(seq_along(cas_type), \(i)cas_type[[i]] %>%
      dplyr::rename(setNames("n", paste0("multi_", i))))
    cas_plotdat <- Reduce(dplyr::full_join, cas_type)
    cas_plotdat <- data.frame(cas_plotdat[, -1], row.names = cas_plotdat[, 1])
    cas_plotdat[is.na(cas_plotdat)] <- 0
    # p=pcutils::stackplot(cas_plotdat)
    p <- do.call(\(...)pcutils::stackplot(cas_plotdat, ...), plot_param)
    return(p)
  }
  cas_plotdat <- dplyr::group_by(x[, c("type", "subtype", "n")], type, subtype) %>%
    dplyr::summarise(n = sum(n)) %>%
    as.data.frame()

  if (mode == 1) p <- do.call(\(...)pcutils::gghuan2(cas_plotdat, ...), plot_param)
  if (mode == 2) p <- do.call(\(...)pcutils::gghuan2(cas_plotdat[, c("subtype", "type", "n")], ...), plot_param)
  if (mode == 3) {
    p <- do.call(
      \(...)pcutils::my_sankey(cas_plotdat[, c("subtype", "type", "n")], ...),
      pcutils::update_param(list(mode = "gg", num = TRUE), plot_param)
    )
  }
  # if(mode==4)p=pcutils::my_circo(cas_plotdat,...)

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
show_cas_type <- function() {
  pcutils::lib_ps("ggtree", "gggenes", "aplot", library = FALSE)
  # load("cas_info.rda",envir = environment())
  data("cas_info", package = "iCRISPR", envir = environment())
  # cas_type_df$cas1=strsplit(cas_type_df$cas,"[,-]")%>%sapply(.,\(i)i[1])

  cas_tree <- ggtree::fortify(cas_tree)
  tree <- ggtree::ggtree(cas_tree, ladderize = FALSE) +
    ggtree::geom_tiplab() +
    ggtree::geom_highlight(node = 24, fill = "skyblue", alpha = 0.5, to.bottom = TRUE) +
    ggtree::geom_highlight(node = 25, fill = "pink", alpha = 0.5, to.bottom = TRUE) +
    geom_text(
      data = data.frame(
        x = c(1, 1), y = c(21.5, 9.5),
        label = c("Class 1", "Class 2")
      ),
      mapping = aes(x = x, y = y, label = label)
    ) + xlim(0, 5)

  p <- ggplot(cas_type_df, aes(y = sub_type, fill = cas1)) +
    gggenes::geom_gene_arrow(aes(xmin = position, xmax = position + 1, forward = (strand == "+"))) +
    gggenes::geom_gene_label(aes(xmin = position, xmax = position + 1, label = cas)) +
    pcutils::scale_fill_pc(n = 31)+
    geom_hline(yintercept = 10.5, linetype = 2) +
    theme_void() +
    theme(legend.position = "none")
  aplot::insert_left(p, tree, width = 0.3)
  # ggsave("inst/extdata/show_cas_type.pdf",width = 14,height = 6)
}

#' Summary crispr and spacer number at different evidence levels
#'
#' @param crispr crispr object
#' @param each_genome summary each genome information in multi_crispr?
#' @param use_CCF use_CCF=TRUE means use the spacer number count by CCF, but sometimes do'not match spacer.
#' @return evidence_level object
#' @export
#'
#' @examples
#' level_res <- summary_levels(multi_crispr, each_genome = FALSE)
#' plot(level_res)
summary_levels <- function(crispr, each_genome = FALSE, use_CCF = TRUE) {
  if (inherits(crispr, "multi_crispr")) {
    aaa <- lapply(crispr, summary_levels)
    bbb <- do.call(rbind, aaa)
    rownames(bbb) <- NULL
    if (!each_genome) {
      ccc <- dplyr::group_by(bbb, evidence_level, Cas) %>%
        dplyr::summarise(genome_num = dplyr::n(), crispr_num = sum(crispr_num), spacer_num = sum(spacer_num)) %>%
        as.data.frame()
      ccc <- dplyr::mutate(ccc, spacer_num_per_array = spacer_num / crispr_num)
    } else {
      ccc <- bbb
    }
    class(ccc) <- c("evidence_level", class(ccc))
    return(ccc)
  }
  if (!inherits(crispr, "crispr")) {
    return(NULL)
  }
  if (is.null(crispr$Spacer)) {
    return(NULL)
  }

  if ("evidence_Level" %in% colnames(crispr$CRISPR)) crispr$CRISPR <- dplyr::rename(crispr$CRISPR, evidence_level = "evidence_Level")
  # 有一个bug是CCF结果的spacer_number和真正的GFF给出的spacer数量不一致（罕见错误）
  # 见crisprfinder/output/dir_296/GCA_005025685.1_PDT000277779.2_genomic.fna
  # 见crisprfinder//output/dir_345/GCA_021349275.1_PDT001215381.1_genomic.fna
  if ("spacer_number" %in% colnames(crispr$CRISPR) & use_CCF) {
    n_spacer2 <- dplyr::select(crispr$CRISPR, genome, evidence_level, spacer_number) %>%
      dplyr::group_by(genome, evidence_level) %>%
      dplyr::summarise(crispr_num = dplyr::n(), spacer_num = sum(spacer_number)) %>%
      as.data.frame()
  } else {
    # 旧版本不包含spacer_number，需要统计
    # 发现spacer数量不一致，应该以gff为准，因为我们能做的就是提出来的这部分
    # 所以当use_CCF=FALSE,会慢一些，但是是准确的
    n_spacer <- dplyr::count(crispr$Spacer, genome, CRISPR_id)
    n_spacer2 <- dplyr::select(crispr$CRISPR, genome, CRISPR_id, evidence_level) %>%
      dplyr::left_join(., n_spacer, by = c("genome" = "genome", "CRISPR_id" = "CRISPR_id")) %>%
      dplyr::group_by(genome, evidence_level) %>%
      dplyr::summarise(crispr_num = dplyr::n(), spacer_num = sum(n)) %>%
      as.data.frame()
  }
  n_spacer2$Cas <- ifelse(is.null(crispr$Cas), "Non-Cas", "With-Cas")

  class(n_spacer2) <- c("evidence_level", class(n_spacer2))
  return(n_spacer2)
}

#' Combine some evidence_level object
#'
#' @param ... evidence_level object
#'
#' @return evidence_level object
#' @exportS3Method
#' @method combine evidence_level
#'
combine.evidence_level <- function(...) {
  all_c <- list(...)
  if (!all(sapply(all_c, inherits, what = "evidence_level"))) stop("all input should be evidence_level object")
  res <- lapply(all_c, \(i){
    if ("genome" %in% colnames(i)) {
      ccc <- dplyr::group_by(i, evidence_level, Cas) %>%
        dplyr::summarise(genome_num = dplyr::n(), crispr_num = sum(crispr_num), spacer_num = sum(spacer_num)) %>%
        as.data.frame()
      return(ccc)
    } else {
      return(i)
    }
  })
  all_res <- do.call(rbind, res)
  ccc <- dplyr::group_by(all_res, evidence_level, Cas) %>%
    dplyr::summarise(genome_num = sum(genome_num), crispr_num = sum(crispr_num), spacer_num = sum(spacer_num)) %>%
    as.data.frame()
  ccc <- dplyr::mutate(ccc, spacer_num_per_array = spacer_num / crispr_num)
  class(ccc) <- c("evidence_level", "data.frame")
  ccc
}

#' Plot the distribution of array and spacer number at different evidence levels.
#' @method plot evidence_level
#'
#' @param x evidence_level object
#' @param mode 1~2
#' @param ... additional
#' @param num_size number font size
#'
#' @return ggplot
#' @exportS3Method
#' @rdname summary_levels
plot.evidence_level <- function(x, mode = 1, num_size = 4, ...) {
  pcutils::lib_ps("reshape2", "scales", library = FALSE)
  ccc <- x
  if (TRUE) {
    ccc <- dplyr::group_by(x, evidence_level, Cas) %>%
      dplyr::summarise(crispr_num = sum(crispr_num), spacer_num = sum(spacer_num)) %>%
      as.data.frame()
    ccc <- dplyr::mutate(ccc, spacer_num_per_array = spacer_num / crispr_num)
  }
  ccc$spacer_num_per_array <- round(ccc$spacer_num_per_array, 1)
  plotdat <- reshape2::melt(ccc, id.vars = 1:2)
  if (mode == 1) {
    p <- ggplot(data = plotdat, aes(x = evidence_level, y = value, fill = Cas)) +
      geom_col(position = position_dodge(width = 1)) +
      geom_text(aes(label = value), position = position_dodge(width = 1), vjust = 0, size = num_size) +
      scale_fill_manual(values = c("#a6cee3", "#78c679")) +
      labs(y = "Number") +
      facet_wrap(. ~ variable, nrow = 1, scales = "free_y")
  }
  if (mode == 2) {
    p <- ggplot(data = plotdat, aes(x = evidence_level, y = value, fill = Cas)) +
      geom_bar(stat = "identity", position = position_fill()) +
      scale_y_continuous(labels = scales::percent) +
      geom_text(aes(label = value), position = position_fill(), size = num_size) +
      labs(y = "Percentage") +
      scale_fill_manual(values = c("#a6cee3", "#78c679")) +
      facet_wrap(. ~ variable, nrow = 1)
  }
  p
}
