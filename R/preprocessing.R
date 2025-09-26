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