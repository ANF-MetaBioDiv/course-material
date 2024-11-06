#' Plot one quality profile per fastq file
#'
#' @return a pdf file containing one quality profile per file
#' @param filefw R1 file names
#' @param filefw R2 file names
#' @param outfile Output pdf file name
#' @export

qualityprofile <- function(filefw, filerev, outfile) {
  pdf(outfile, onefile = TRUE)
  for (i in seq_along(filefw)) {
    f <- dada2::plotQualityProfile(filefw[[i]])
    r <- dada2::plotQualityProfile(filerev[[i]])
    print(f)
    print(r)
  }
  dev.off()
}

#' Primer trimming using cutadapt
#'
#' @return pair of gunzipped fastq files (R1 and R2) of primer trimmed sequences
#' @param forward_files vector of names of forward reads files to trim
#' @param reverse_files vector of names of reverse reads files to trim
#' @param primer_fwd forward primer
#' @param primer_rev reverse primer
#' @param output_dir folder name where to write trimmed sequences
#' @param min_size sequences shorter than min_size after trimming are
#' filtered out
#' @details Input file names have to end with .fastq.gz and contain _R1 or _R2.
#'
#' Cutadapt need to be installed on your system.
#' @export

primer_trim <- function(
  forward_files, reverse_files,
  primer_fwd, primer_rev,
  output_dir, min_size
) {

  a <- sub("_R[12].+fastq.gz", "", forward_files)
  b <- sub("_R[12].+fastq.gz", "", reverse_files)

  if (!identical(a, b)) {
    stop("forward_files and reverse_files are not matching")
  }

  path_to_dada2_log <- here::here(
    "outputs",
    "dada2",
    "log"
  )

  if (!dir.exists(path_to_dada2_log)) dir.create(path_to_dada2_log)

  invisible(
    mapply(
      function(x, y) {
        run_cutadapt_command(
          forward_file = x,
          reverse_file = y,
          primer_fwd = primer_fwd,
          primer_rev = primer_rev,
          output_dir = output_dir,
          min_size = min_size,
          log_dir = path_to_dada2_log
        )
      },
      forward_files,
      reverse_files
    )
  )

  all_logs_names <- list.files(
    path_to_dada2_log,
    pattern = "_primer_trimming.log",
    full.names = TRUE
  )

  all_logs <- lapply(all_logs_names, read.table, sep = "\t")

  all_logs <- Reduce("rbind", all_logs)

  names(all_logs) <- c(
    "sample", "status", "in_reads", "in_bp",
    "too_short", "too_long", "too_many_n",
    "out_reads", "w/adapters", "qualtrim_bp",
    "out_bp", "w/adapters2", "qualtrim2_bp", "out2_bp"
  )

  log_file_export <- file.path(path_to_dada2_log, "primer_trimming.log")

  write.table(all_logs, log_file_export, row.names = FALSE, quote = FALSE)

  invisible(file.remove(all_logs_names))

  return(all_logs)
}

#' Run cutadapt (Unix)
#' @description R function to call cutadapt on the system (Unix) to trim one sample
#' @keywords internal

run_cutadapt_command <- function(
  forward_file, reverse_file,
  primer_fwd, primer_rev,
  output_dir, min_size, log_dir
) {

  bash_script <- here::here("bash", "primer_trimming.bash")

  command <- paste(
    "bash",
    bash_script,
    forward_file,
    reverse_file,
    primer_fwd,
    primer_rev,
    output_dir,
    min_size,
    log_dir
  )
  system(command)
}

#' Primer trimming using cutadapt on windows
#'
#' @return pair of gunzipped fastq files (R1 and R2) of primer trimmed sequences
#' @param forward_files vector of names of forward reads files to trim
#' @param reverse_files vector of names of reverse reads files to trim
#' @param primer_fwd forward primer
#' @param primer_rev reverse primer
#' @param output_dir folder name where to write trimmed sequences
#' @param min_size sequences shorter than min_size after trimming are
#' filtered out
#' @details Input file names have to end with .fastq.gz and contain _R1 or _R2
#'
#' Cutadapt has to be installed on your system.
#' @export

primer_trim_windows <- function(
  forward_files, reverse_files,
  primer_fwd, primer_rev,
  output_dir, min_size
) {

  a <- sub("_R[12].+fastq.gz", "", forward_files)
  b <- sub("_R[12].+fastq.gz", "", reverse_files)

  if (!identical(a, b)) {
    stop("forward_files and reverse_files are not matching")
  }

  path_to_dada2_log <- here::here(
    "outputs",
    "dada2",
    "log"
  )

  if (!dir.exists(path_to_dada2_log)) dir.create(path_to_dada2_log)

  invisible(
    mapply(
      function(x, y) {
        run_cutadapt_command_windows(
          forward_file = x,
          reverse_file = y,
          primer_fwd = primer_fwd,
          primer_rev = primer_rev,
          output_dir = output_dir,
          min_size = min_size,
          log_dir = path_to_dada2_log,
          sample_name = basename(sub("_R[12].+fastq.gz", "", x))
        )
      },
      forward_files,
      reverse_files
    )
  )

  all_logs_names <- list.files(
    path_to_dada2_log,
    pattern = "_primer_trimming.log",
    full.names = TRUE
  )

  all_logs <- lapply(all_logs_names, read.table, sep = "\t")

  all_logs <- Reduce("rbind", all_logs)

  names(all_logs) <- c(
    "sample", "status", "in_reads", "in_bp",
    "too_short", "too_long", "too_many_n",
    "out_reads", "w/adapters", "qualtrim_bp",
    "out_bp", "w/adapters2", "qualtrim2_bp", "out2_bp"
  )

  log_file_export <- file.path(path_to_dada2_log, "primer_trimming.log")

  write.table(all_logs, log_file_export, row.names = FALSE, quote = FALSE)

  invisible(file.remove(all_logs_names))

  return(all_logs)
}

#' Run cutadapt (Windows)
#' @description R function to call cutadapt on the system (Windows) to trim one sample
#' @keywords internal

run_cutadapt_command_windows <- function(
  forward_file, reverse_file,
  primer_fwd, primer_rev,
  output_dir, min_size, log_dir, sample_name
) {

  command <- paste(
    "py -m cutadapt",
    "-g", primer_fwd,
    "-G", primer_rev,
    "--report=minimal",
    "--discard-untrimmed",
    "--minimum-length", min_size,
    "--no-indels",
    "-o", here::here(output_dir, basename(forward_file)),
    "-p", here::here(output_dir, basename(reverse_file)),
    forward_file,
    reverse_file
  )

  log <- system(command, intern = TRUE)

  log <- paste(sample_name, log[2], sep = "\t")

  cat(
    log,
    file = here::here(log_dir, paste0(sample_name, "_primer_trimming.log")),
    sep = "\n"
  )

}