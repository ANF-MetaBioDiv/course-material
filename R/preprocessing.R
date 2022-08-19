
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