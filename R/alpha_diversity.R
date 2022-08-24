#' Plot rarefaction curves
#'
#' The function plot rarefaction curves...
#'
#' @return Print ggplot rarefaction curves
#' @param physeq A phyloseq class object, from which abundance
#' data are extracted
#' @param step Step size for sample size in rarefaction curves
#' @param label Default `NULL`. Character string.
#' The name of the variable to map to text labels on the plot.
#' Similar to color option but for plotting text.
#' @param color (Optional). Default ‘NULL’. Character string.
#' The name of the variable to map to colors in the plot.
#' This can be a sample variable (among the set returned by
#' ‘sample_variables(physeq)’ ) or taxonomic rank (among the set returned by
#' ‘rank_names(physeq)’).
#' The color scheme is chosen automatically by
#' ‘link{ggplot}’, but it can be modified afterward with an
#' additional layer using ‘scale_color_manual’.
#' @param plot Logical, should the graphic be plotted.
#' @param parallel should rarefaction be parallelized (using parallel framework)
#' @param se Default TRUE. Logical. Should standard errors be computed.

ggrare <- function(physeq, step = 10, label = NULL, color = NULL,
                   plot = TRUE, parallel = FALSE, se = TRUE) {

  require("ggplot2")

  x <- as(phyloseq::otu_table(physeq), "matrix")
  if (phyloseq::taxa_are_rows(physeq)) x <- t(x)
  ## This script is adapted from vegan `rarecurve` function
  tot <- rowSums(x)
  S <- rowSums(x > 0)
  nr <- nrow(x)

  rarefun <- function(i) {
    cat(paste("rarefying sample", rownames(x)[i]), sep = "\n")
    n <- seq(1, tot[i], by = step)
    if (n[length(n)] != tot[i]) {
      n <- c(n, tot[i])
    }
    y <- vegan::rarefy(x[i, ,drop = FALSE], n, se = se)
    if (nrow(y) != 1) {
      rownames(y) <- c(".S", ".se")
      return(data.frame(t(y), Size = n, Sample = rownames(x)[i]))
    } else {
      return(data.frame(.S = y[1, ], Size = n, Sample = rownames(x)[i]))
    }
  }

  if (parallel) {
    out <- parallel::mclapply(seq_len(nr), rarefun, mc.preschedule = FALSE)
  } else {
    out <- lapply(seq_len(nr), rarefun)
  }

  df <- do.call(rbind, out)

  ## Get sample data 
  if (!is.null(phyloseq::sample_data(physeq, FALSE))) {
    sdf <- as(phyloseq::sample_data(physeq), "data.frame")
    sdf$Sample <- rownames(sdf)
    data <- merge(df, sdf, by = "Sample")
    labels <- data.frame(x = tot, y = S, Sample = rownames(x))
    labels <- merge(labels, sdf, by = "Sample")
  }

  ## Add, any custom-supplied plot-mapped variables
  if (length(color) > 1) {
    data$color <- color
    names(data)[names(data) == "color"] <- deparse(substitute(color))
    color <- deparse(substitute(color))
  }

  if (length(label) > 1) {
    labels$label <- label
    names(labels)[names(labels) == "label"] <- deparse(substitute(label))
    label <- deparse(substitute(label))
  }

  p <- ggplot(data = data,
              aes_string(x = "Size", y = ".S",
                         group = "Sample", color = color)) +
    labs(x = "Sample Size", y = "Species Richness")

  if (!is.null(label)) {
    p <- p + geom_text(data = labels,
                      aes_string(x = "x", y = "y",
                                 label = label, color = color),
                      size = 4, hjust = 0)
  }

  p <- p + geom_line()

  if (se) { ## add standard error if available
    p <- p +
        geom_ribbon(aes_string(ymin = ".S - .se", ymax = ".S + .se",
                               color = NULL, fill = color),
                    alpha = 0.2)
  }

  if (plot) {
    plot(p)
  }

  invisible(p)

}

#' Check for normal distribution of alpha diversity measurements
#'
#'
#' @return Print ggplot rarefaction curves
#' @param rich A dataframe with alpha diversity indices as columns
#' and samples as rows as the output of pyloseq::estimate_richness()
#' @param file Output file name

indices_normality <- function(rich, nrow, ncol) {

  ### p-value < 0.05 means data failed normality test

  par(mfrow = c(nrow, ncol))

  for (i in names(rich)) {
    shap <- shapiro.test(rich[, i])
    qqnorm(rich[, i], main = i, sub = shap$p.value)
    qqline(rich[, i])
  }

  par(mfrow = c(1, 1))
}