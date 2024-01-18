install_blast <- function(blast_tar_gz = NULL, blast_version = "ncbi-blast-2.14.0+") {
    # Detect the operating system
    os <- tolower(Sys.info()[["sysname"]])
    # Set the installation URL and command based on the operating system
    if (os == "linux") {
        url <- paste0("https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/", blast_version, "-x64-linux.tar.gz")
        install_cmd <- "tar -xf ncbi-blast*.tar.gz"
    } else if (os == "darwin") {
        url <- paste0("https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/", blast_version, "-x64-macosx.tar.gz")
        install_cmd <- "tar -xf ncbi-blast*.tar.gz"
    } else if (os == "windows") {
        url <- paste0("https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/", blast_version, "-x64-win64.tar.gz")
        install_cmd <- "unzip ncbi-blast*.tar.gz \n ren ncbi-blast* ncbi-blast"
    } else {
        stop("Unsupported operating system!")
    }
    # Set the file name and installation directory
    filename <- basename(url)

    install_dir <- file.path(system.file(package = "iCRISPR"), "software")
    # Create the installation directory if it does not exist
    dir.create(install_dir, recursive = TRUE, showWarnings = FALSE)

    # Set the full path to the installation directory
    install_path <- file.path(install_dir, filename)

    if (is.null(blast_tar_gz)) {
        # Download the file
        ori_time <- getOption("timeout")
        options(timeout = 300)
        tryCatch(expr = {
            download.file(url, destfile = install_path)
        }, error = function(e) {
            options(timeout = ori_time)
            stop("Try download yourself from https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/")
        })
    } else {
        if (file.exists(blast_tar_gz) & grepl("ncbi-blast.*\\.tar\\.gz", blast_tar_gz)) {
            filename <- basename(blast_tar_gz)
            install_path <- file.path(install_dir, filename)
            file.copy(blast_tar_gz, install_path)
        } else {
            stop("Wrong file: ", blast_tar_gz)
        }
    }

    # Extract the downloaded file
    setwd(install_dir) # Change to the installation directory
    system(install_cmd)

    # Remove the downloaded file
    file.remove(install_path)

    message("BLAST has been successfully installed!\n")
}

run_blast <- function(query_file, database_file, output_file, blast_program = "blastn", num_threads = 1) {
    # 构建BLAST命令
    blast_cmd <- paste0(blast_program, " -query ", query_file, " -db ", database_file, " -out ", output_file, " -num_threads ", num_threads)

    # 执行BLAST命令
    system(blast_cmd)
}

# ==========blast 后用taxonkit lca 确定物种，或者注释功能

#' Generate Random Sequences
#'
#' @param num_seqs The number of sequences to generate.
#' @param mean_length The mean length of the sequences.
#' @param sd_length The standard deviation of the sequence lengths.
#' @param element default: c("A", "T", "C", "G")
#'
#' @return A data frame with sequence names and sequences.
#' @export
#' @examples
#' random_seq()
#'
random_seq <- function(num_seqs = 100, mean_length = 32, sd_length = 3, element = c("A", "T", "C", "G")) {
    DNA <- c("A", "T", "C", "G")
    Protein <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")

    seq_length <- abs(round(rnorm(num_seqs, mean = mean_length, sd = sd_length)))
    seq_length[seq_length < 1] <- 1

    seqs <- lapply(seq_len(num_seqs), \(i){
        seq <- paste(sample(element, seq_length[i], replace = TRUE), collapse = "")
    })
    sequences <- data.frame(
        Sequence_ID = paste0("random", seq_len(num_seqs)),
        Sequence = unlist(seqs),
        stringsAsFactors = FALSE
    )
    return(sequences)
}

#' Generate Random Tandem Repeat Sequences
#'
#' @param num_seqs The number of sequences to generate.
#' @param mean_unit_length The mean length of the repeat unit.
#' @param sd_unit_length The standard deviation of the repeat unit lengths.
#' @param mean_num_units The mean number of repeat units in each sequence.
#' @param sd_num_units The standard deviation of the number of repeat units.
#' @param element Default: c("A", "T", "C", "G")
#'
#' @return A data frame with sequence names and tandem repeat sequences.
#' @export
#' @examples
#' random_tandem_repeat()
random_tandem_repeat <- function(num_seqs = 100, mean_unit_length = 6, sd_unit_length = 1,
                                 mean_num_units = 5, sd_num_units = 2, element = c("A", "T", "C", "G")) {
    # Generate random repeat unit using random_seq function
    repeat_unit <- random_seq(num_seqs,
        mean_length = mean_unit_length,
        sd_length = sd_unit_length, element = element
    )$Sequence
    # Generate random number of repeat units
    num_units <- abs(round(rnorm(num_seqs, mean = mean_num_units, sd = sd_num_units), 2))
    num_units[num_units < 1] <- 1

    seqs <- vapply(seq_len(num_seqs), \(i){
        seq <- paste(rep(repeat_unit[i], ceiling(num_units[i])), collapse = "")
        len <- round(num_units[i] * nchar(repeat_unit[i]))
        seq <- substr(seq, 1, len)
        return(seq)
    }, FUN.VALUE = character(1))

    sequences <- data.frame(
        seq_name = paste0("tandem_repeat", seq_len(num_seqs)),
        seq = seqs,
        repeat_unit = repeat_unit,
        num_units = num_units,
        stringsAsFactors = FALSE
    )

    return(sequences)
}

#' Count kmers in a sequence
#'
#' @param sequence The input DNA sequence.
#' @param k_size The length of the kmers.
#' @export
#' @return A named vector with kmer counts.
count_kmers <- function(sequence, k_size) {
    data <- numeric()
    size <- nchar(sequence)

    for (i in 1:(size - k_size + 1)) {
        kmer <- substr(sequence, i, i + k_size - 1)
        if (is.element(kmer, names(data))) {
            data[kmer] <- data[kmer] + 1
        } else {
            data[kmer] <- 1
        }
    }

    return(data)
}

#' Calculate Sequence Complexity of a sequence
#'
#' @param sequence The input DNA sequence.
#' @param max_k The maximum k value for kmers.
#' @return The Sequence Complexity value.
#' @export
#' @references https://resources.qiagenbioinformatics.com/manuals/clccancerresearchworkbench/200/index.php?manual=How_sequence_complexity_is_calculated.html
cal_sc <- function(sequence, max_k = NULL) {
    fenmu <- function(i, lens) {
        a <- 4^i
        b <- lens + 1 - i
        if (a > b) {
            a <- b
        }
        return(a)
    }
    sc <- 1
    L <- nchar(sequence)
    if (is.null(max_k)) {
        max_k <- 7
    }
    max_k <- min(max_k, L)

    for (k in 1:(max_k - 1)) {
        sc <- sc * length(count_kmers(sequence, k)) / fenmu(k, L)
    }
    return(round(sc, 5))
}

#' Calculate Shannon Complexity of a sequence
#'
#' @param sequence The input DNA sequence.
#' @param base default the number of elements type in sequence.
#' @export
#' @return The Shannon Complexity value.
#' @references 1. Konopka, A. K. Sequence Complexity and Composition. in eLS (ed. John Wiley & Sons, Ltd) (Wiley, 2005). doi:10.1038/npg.els.0005260.
cal_shannon <- function(sequence, base = NULL) {
    x <- count_kmers(sequence, 1)
    if (is.null(base)) base <- length(x)
    x <- x / sum(x)
    sum(-x * log(x, base))
}

#' Summarize sequence information
#'
#' @param sequence The input DNA sequence.
#' @param seq_col the column name of sequence if your input is a dataframe.
#' @param max_k The maximum k value for kmers.
#' @param threads threads
#' @param verbose verbose
#'
#' @return data.frame
#' @export
#' @examples
#' summary_seq("CCTGAACCTATGCCGTCCACCTTGCGTTGCCTT")
#' summary_seq(random_seq(10))
summary_seq <- function(sequence, max_k = 7, seq_col = 2, threads = 1, verbose = TRUE) {
    summary_seq_in <- function(sequence, max_k = 7) {
        # devtools::load_all("~/Documents/R/iCRISPR/")
        counts <- count_kmers(sequence, 1) %>% as.list()
        for (i in c("A", "T", "C", "G")) {
            if (is.null(counts[[i]])) counts[[i]] <- 0
        }
        counts$shannon <- round(cal_shannon(sequence), 4)
        counts$complexity <- round(cal_sc(sequence, max_k), 4)
        res <- data.frame("A" = counts$A, "T" = counts$`T`, "C" = counts$C, "G" = counts$G, "shannon" = counts$shannon, "complexity" = counts$complexity)
        res[is.na(res)] <- 0
        res <- dplyr::mutate(res, length = nchar(sequence), GC_content = round((G + C) / length, 4))
        return(res)
    }
    flag <- is.data.frame(sequence)
    if (flag) {
        df <- sequence
        sequence <- sequence[, seq_col, drop = TRUE]
    }

    # parallel
    reps <- length(sequence)
    # main function
    loop <- function(i) {
        summary_seq_in(sequence[i], max_k = max_k)
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
            res <- lapply(1:reps, loop)
        }
    }
    # simplify method
    res <- do.call(rbind, res)

    if (flag) {
        rownames(res) <- NULL
        return(cbind(df, res))
    } else {
        return(cbind(sequence = sequence, res))
    }
}
