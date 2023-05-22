install_blast <- function(blast_tar_gz=NULL,blast_version = "ncbi-blast-2.14.0+") {
  # Detect the operating system
  os <- tolower(Sys.info()[["sysname"]])
  # Set the installation URL and command based on the operating system
  if (os == "linux") {
    url <- paste0("https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/",blast_version,"-x64-linux.tar.gz")
    install_cmd <- "tar -xf ncbi-blast*.tar.gz"
  } else if (os == "darwin") {
    url <- paste0("https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/",blast_version,"-x64-macosx.tar.gz")
    install_cmd <- "tar -xf ncbi-blast*.tar.gz"
  } else if (os == "windows") {
    url <- paste0("https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/",blast_version,"-x64-win64.tar.gz")
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

  if(is.null(blast_tar_gz)){
    # Download the file
    ori_time=getOption("timeout")
    options(timeout = 300)
    tryCatch(expr = {
      download.file(url, destfile = install_path)
    },error=function(e){
      options(timeout = ori_time);
      stop("Try download yourself from https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/")
    })
  }
  else {
    if(file.exists(blast_tar_gz)&grepl("ncbi-blast.*\\.tar\\.gz",blast_tar_gz)){
      filename <- basename(blast_tar_gz)
      install_path <- file.path(install_dir, filename)
      file.copy(blast_tar_gz,install_path)
    }
    else stop("Wrong file: ",blast_tar_gz)
  }

  # Extract the downloaded file
  setwd(install_dir)  # Change to the installation directory
  system(install_cmd)

  # Remove the downloaded file
  file.remove(install_path)

  cat("BLAST has been successfully installed!\n")
}

run_blast <- function(query_file, database_file, output_file, blast_program = "blastn", num_threads = 1) {
  # 构建BLAST命令
  blast_cmd <- paste0(blast_program, " -query ", query_file, " -db ", database_file, " -out ", output_file, " -num_threads ", num_threads)

  # 执行BLAST命令
  system(blast_cmd)
}
