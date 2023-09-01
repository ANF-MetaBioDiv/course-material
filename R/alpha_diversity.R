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

rankabuncomp<-function (x, y = NULL, factor = NULL, return.data = T, specnames = c(1:3), 
                        scale = "abundance", scaledx = F, type = "o", rainbow = T, 
                        legendpos = "topright", xlim = c(1, max1), ylim = c(0, max2), ...) 
{
  groups <- table(y[, factor])
  levels <- names(groups)
  m <- length(groups)
  max1 <- max(diversitycomp(x, y, factor1 = factor, index = "richness", 
                            method = "pooled")[, 2])
  if (scaledx == T) {
    xlim <- c(0, 100)
  }
  max2 <- max.2 <- 0
  for (i in 1:m) {
    if (scale == "abundance") {
      max.2 <- rankabundance(x, y, factor, levels[i])[1, 
                                                      "abundance"]
    }
    if (scale == "logabun") {
      max.2 <- rankabundance(x, y, factor, levels[i])[1, 
                                                      "abundance"]
    }
    if (scale == "proportion") {
      max.2 <- rankabundance(x, y, factor, levels[i])[1, 
                                                      "proportion"]
    }
    if (max.2 > max2) {
      max2 <- max.2
    }
  }
  if (scale == "accumfreq") {
    max2 <- 100
  }
  max2 <- as.numeric(max2)
  if (rainbow == F) {
    if (scale == "logabun" && all.equal(ylim, c(0, max2)) == 
        T) {
      ylim <- c(1, max2)
    }
    rankabunplot(rankabundance(x, y, factor, levels[1]), 
                 scale = scale, scaledx = scaledx, type = type, labels = levels[1], 
                 xlim = xlim, ylim = ylim, pch = 1, specnames = NULL, 
                 ...)
    for (i in 2:m) {
      rankabunplot(rankabundance(x, y, factor, levels[i]), 
                   addit = T, scale = scale, scaledx = scaledx, 
                   type = type, labels = levels[i], pch = i, specnames = NULL, 
                   ...)
    }
    legend(legendpos, legend = levels, pch = c(1:m))
  }
  else {
    grDevices::palette(colorspace::rainbow_hcl(m, c = 90, 
                                               l = 50))
    if (scale == "logabun" && all.equal(ylim, c(0, max2)) == 
        T) {
      ylim <- c(1, max2)
    }
    rankabunplot(rankabundance(x, y, factor, levels[1]), 
                 scale = scale, scaledx = scaledx, type = type, labels = levels[1], 
                 xlim = xlim, ylim = ylim, col = 1, pch = 1, specnames = NULL, 
                 ...)
    for (i in 2:m) {
      rankabunplot(rankabundance(x, y, factor, levels[i]), 
                   addit = T, scale = scale, scaledx = scaledx, 
                   type = type, labels = levels[i], col = i, pch = i, 
                   specnames = NULL, ...)
    }
    legend(legendpos, legend = levels, pch = c(1:m), 
           col = c(1:m))
    grDevices::palette("default")
  }
  if (return.data == T) {
    for (i in 1:m) {
      resulti <- data.frame(rankabundance(x, y, factor, 
                                          levels[i]))
      resulti <- data.frame(Grouping = rep(levels[i], nrow(resulti)), 
                            species = rownames(resulti), labelit = rep(FALSE, 
                                                                       nrow(resulti)), resulti)
      spec.max <- min(max(specnames), nrow(resulti))
      resulti[c(1:spec.max), "labelit"] <- as.logical(1)
      rownames(resulti) <- NULL
      if (i == 1) {
        result <- resulti
      }
      else {
        result <- rbind(result, resulti)
      }
    }
    return(result)
  }
}


rankabunplot<-function (xr, addit = F, labels = "", scale = "abundance", scaledx = F, 
                        type = "o", xlim = c(min(xpos), max(xpos)), ylim = c(0, max(x[, 
                                                                                      scale])), specnames = c(1:5), srt = 0, ...) 
{
  x <- xr
  xpos <- 1:nrow(x)
  if (scaledx == T) {
    xpos <- xpos/nrow(x) * 100
  }
  if (scale == "accumfreq") {
    type <- "o"
  }
  if (addit == F) {
    if (scale == "logabun") {
      if (all.equal(ylim, c(0, max(x[, scale]))) == T) {
        ylim <- c(1, max(x[, "abundance"]))
      }
      graphics::plot(xpos, x[, "abundance"], xlab = "species rank", 
                     ylab = "abundance", type = type, bty = "l", log = "y", 
                     xlim = xlim, ylim = ylim, ...)
    }
    else {
      graphics::plot(xpos, x[, scale], xlab = "species rank", 
                     ylab = scale, type = type, bty = "l", ylim = ylim, 
                     xlim = xlim, ...)
    }
  }
  else {
    if (scale == "logabun") {
      graphics::points(xpos, x[, "abundance"], type = type, 
                       ...)
    }
    else {
      graphics::points(xpos, x[, scale], type = type, ...)
    }
  }
  if (length(specnames) > 0) {
    names.space <- paste0("  ", rownames(x))
    for (i in specnames) {
      if (scale == "logabun") {
        graphics::text(i, x[i, "abundance"], names.space[i], 
                       pos = 4, srt = srt, offset = 0, adj = 1)
      }
      else {
        graphics::text(i, x[i, scale], names.space[i], 
                       pos = 4, srt = srt, offset = 0, adj = 1)
      }
    }
  }
  if (labels != "") {
    if (scale == "logabun") {
      graphics::text(1, x[1, "abundance"], labels = labels, 
                     pos = 2)
    }
    else {
      graphics::text(1, x[1, scale], labels = labels, pos = 2)
    }
  }
}

rankabundance<-function (x, y = "", factor = "", level, digits = 1, t = qt(0.975, 
                                                                           df = n - 1)) 
{
  if (inherits(y, "data.frame") && factor != "") {
    subs <- y[, factor] == level
    for (q in 1:length(subs)) {
      if (is.na(subs[q])) {
        subs[q] <- F
      }
    }
    x <- x[subs, , drop = F]
    freq <- apply(x, 2, sum)
    subs <- freq > 0
    x <- x[, subs, drop = F]
  }
  if (dim(as.matrix(x))[1] == 0) {
    result <- array(NA, dim = c(1, 8))
    colnames(result) <- c("rank", "abundance", "proportion", 
                          "plower", "pupper", "accumfreq", "logabun", "rankfreq")
    rownames(result) <- "none"
    return(result)
  }
  total <- apply(x, 1, sum)
  p <- ncol(x)
  n <- nrow(x)
  mu <- sum(total)/n
  result <- array(dim = c(p, 8))
  colnames(result) <- c("rank", "abundance", "proportion", 
                        "plower", "pupper", "accumfreq", "logabun", "rankfreq")
  rownames(result) <- colnames(x)
  for (j in 1:p) {
    spec <- x[, j]
    pi <- spec/total
    p <- sum(spec)/sum(total)
    sigma2 <- 0
    for (i in 1:n) {
      sigma2 <- sigma2 + (total[i]^2 * (pi[i] - p)^2)
    }
    sigma2 <- sigma2/(n * (n - 1) * mu * mu)
    sigma <- sigma2^0.5
    result[j, 2] <- sum(spec)
    result[j, 3] <- p * 100
    result[j, 4] <- (p - t * sigma) * 100
    result[j, 5] <- (p + t * sigma) * 100
  }
  p <- ncol(x)
  result2 <- result
  seq <- rev(order(result[, 2], -order(rownames(result))))
  result[1:p, ] <- result2[seq, ]
  rownames(result)[1:p] <- rownames(result2)[seq]
  result[, 1] <- c(1:ncol(x))
  result[, 6] <- cumsum(result[, 3])
  result[, 7] <- log(result[, 2], base = 10)
  result[, 8] <- result[, 1]/ncol(x) * 100
  result <- round(result, digits = digits)
  return(result)
}

diversitycomp<-function (x, y = NULL, factor1 = NULL, factor2 = NULL, index = c("Shannon", 
                                                                                "Simpson", "inverseSimpson", "Logalpha", "Berger", "simpson.unb", 
                                                                                "simpson.unb.inverse", "richness", "abundance", "Jevenness", 
                                                                                "Eevenness", "jack1", "jack2", "chao", "boot"), method = c("pooled", 
                                                                                                                                           "mean", "sd", "max", "jackknife"), sortit = FALSE, digits = 8) 
{
  INDEX <- c("Shannon", "Simpson", "inverseSimpson", "Logalpha", 
             "Berger", "simpson.unb", "simpson.unb.inverse", "richness", 
             "abundance", "Jevenness", "Eevenness", "jack1", "jack2", 
             "chao", "boot")
  if ((index %in% INDEX) == F) {
    stop(paste("choose an accepted index, not index: ", index, 
               sep = ""))
  }
  METHOD <- c("pooled", "mean", "sd", "max", "jackknife")
  if ((method %in% METHOD) == F) {
    stop(paste("choose an accepted method, not method: ", 
               method, sep = ""))
  }
  if (is.null(y) == F) {
    if ((factor1 %in% names(y)) == F) {
      stop("specified factor1 '", factor1, "' is not a variable of the environmental data frame")
    }
    if (is.factor(y[, factor1]) == F) {
      stop("specified factor1 '", factor1, "' is not a factor")
    }
    y[, factor1] <- as.factor(as.character(y[, factor1]))
    if (is.null(factor2) == F) {
      if ((factor2 %in% names(y)) == F) {
        stop("specified factor2 '", factor2, "' is not a variable of the environmental data frame")
      }
      if (is.factor(y[, factor2]) == F) {
        stop("specified factor2 '", factor2, "' is not a factor")
      }
      y[, factor2] <- as.factor(as.character(y[, factor2]))
    }
  }
  if (is.null(factor2) == T) {
    groups <- table(y[, factor1])
    m <- length(groups)
    levels <- as.character(names(groups))
    result <- array(NA, dim = c(m, 2))
    result[, 1] <- groups
    dimnames(result) <- list(factor1 = levels, c("n", index))
    names(dimnames(result)) <- c(factor1, "")
    for (i in 1:m) {
      if (method %in% c("pooled", "mean", "max", "sd")) {
        result[i, 2] <- as.numeric(diversityresult(x, 
                                                   y, factor = factor1, level = levels[i], method = method, 
                                                   index = index, digits = digits))
      }
      if (method == "jackknife") {
        resultx <- diversityresult(x, y, factor = factor1, 
                                   level = levels[i], method = "jackknife", index = index, 
                                   digits = digits)
        result[i, 2] <- as.numeric(resultx$jack.estimate)
      }
    }
    if (sortit == T) {
      result2 <- result
      seq <- order(result[, 2])
      for (i in 1:m) {
        result[1:m, ] <- result2[seq, ]
      }
      rownames(result) <- rownames(result2)[seq]
    }
    return(result)
  }
  else {
    if (method == "jackknife") {
      stop("jackknife analysis problematic with two factors")
    }
    groups <- table(y[, factor1], y[, factor2])
    levels1 <- rownames(groups)
    levels2 <- colnames(groups)
    m1 <- length(levels1)
    m2 <- length(levels2)
    result <- array(NA, dim = c(m1, m2, 2))
    result[, , 1] <- groups
    dimnames(result) <- list(factor1 = levels1, factor2 = levels2, 
                             c("n", index))
    names(dimnames(result)) <- c(factor1, factor2, "")
    for (i in 1:m1) {
      for (j in 1:m2) {
        if (as.numeric(groups[i, j]) > 0) {
          subs <- y[, factor1] == as.character(levels1[i])
          x1 <- x[subs, , drop = F]
          y1 <- y[subs, , drop = F]
          result[i, j, 2] <- as.numeric(diversityresult(x1, 
                                                        y = y1, factor = factor2, level = levels2[j], 
                                                        method = method, index = index, digits = digits))
        }
      }
    }
    return(result)
  }
}

diversityresult<-function (x, y = NULL, factor = NULL, level = NULL, index = c("Shannon", 
                                                                               "Simpson", "inverseSimpson", "Logalpha", "Berger", "simpson.unb", 
                                                                               "simpson.unb.inverse", "richness", "abundance", "Jevenness", 
                                                                               "Eevenness", "jack1", "jack2", "chao", "boot"), method = c("pooled", 
                                                                                                                                          "each site", "mean", "sd", "max", "jackknife"), sortit = FALSE, 
                           digits = 8) 
{
  INDEX <- c("Shannon", "Simpson", "inverseSimpson", "Logalpha", 
             "Berger", "simpson.unb", "simpson.unb.inverse", "richness", 
             "abundance", "Jevenness", "Eevenness", "jack1", "jack2", 
             "chao", "boot")
  if ((index %in% INDEX) == F) {
    stop(paste("choose an accepted index, not index: ", index, 
               sep = ""))
  }
  METHOD <- c("pooled", "each site", "mean", "sd", "max", "jackknife")
  if ((method %in% METHOD) == F) {
    stop(paste("choose an accepted method, not method: ",
               method, sep = ""))
  }
  if (is.null(y) == F) {
    if ((factor %in% names(y)) == F) {
      stop("specified factor '", factor, "' is not a variable of the environmental data frame")
    }
    if (is.factor(y[, factor]) == F) {
      stop("specified factor '", factor, "' is not a factor")
    }
    levels1 <- as.character(levels(as.factor(as.character(y[, 
                                                            factor]))))
    if ((level %in% levels1) == F) {
      stop("specified level '", level, "' is not an available factor level")
    }
  }
  if (index %in% c("jack1", "jack2", "chao", "boot") && method != "pooled") {
    cat(paste("\n", "Note that default method for index '",
              index, "' is method 'pooled'", "\n\n", sep = ""))
    method <- "pooled"
  }
  diversityresult0 = function(x, index, method) {
    marg <- 1
    if (method == "pooled" && index != "jack1" && index != 
        "jack2" && index != "chao" && index != "boot") {
      x <- apply(x, 2, sum)
      marg <- 2
    }
    if (index == "Shannon") {
      result <- vegan::diversity(x, index = "shannon", MARGIN = marg)
    }
    if (index == "Simpson") {
      result <- vegan::diversity(x, index = "simpson", MARGIN = marg)
    }
    if (index == "inverseSimpson") {
      result <- vegan::diversity(x, index = "invsimpson", MARGIN = marg)
    }
    if (index == "simpson.unb") {
      result <- vegan::simpson.unb(x)
    }
    if (index == "simpson.unb.inverse") {
      result <- vegan::simpson.unb(x, inverse = TRUE)
    }
    if (index == "Logalpha") {
      result <- vegan::fisher.alpha(x, MARGIN = 1)
    }
    if (index == "Berger") {
      if (marg == 2) {
        result <- max(x)/sum(x)
      }
      else {
        tots <- as.matrix(apply(x, marg, sum))
        result <- as.matrix(apply(x, marg, max))
        result <- as.matrix(result/tots)[, 1]
      }
    }
    if (index == "richness") {
      if (marg == 2) {
        result <- sum(x > 0)
      }
      else {
        result <- as.matrix(apply(x > 0, marg, sum))
        result <- result[, 1]
      }
    }
    if (index == "abundance") {
      if (marg == 2) {
        result <- sum(x)
      }
      else {
        result <- as.matrix(apply(x, marg, sum))
        result <- result[, 1]
      }
    }
    if (index == "Jevenness") {
      result1 <- vegan::diversity(x, index = "shannon", MARGIN = marg)
      if (marg == 2) {
        result2 <- sum(x > 0)
      }
      else {
        result2 <- as.matrix(apply(x > 0, marg, sum))
        result2 <- result2[, 1]
      }
      result <- result1/log(result2)
    }
    if (index == "Eevenness") {
      result1 <- vegan::diversity(x, index = "shannon", MARGIN = marg)
      if (marg == 2) {
        result2 <- sum(x > 0)
      }
      else {
        result2 <- as.matrix(apply(x > 0, marg, sum))
        result2 <- result2[, 1]
      }
      result <- exp(result1) / result2
    }
    if (index == "jack1") {
      result <- vegan::specpool(x)$jack1
    }
    if (index == "jack2") {
      result <- vegan::specpool(x)$jack2
    }
    if (index == "chao") {
      result <- vegan::specpool(x)$chao
    }
    if (index == "boot") {
      result <- vegan::specpool(x)$boot
    }
    return(result)
  }
  options(digits = digits)
  if (is.null(y) == F) {
    subs <- y[, factor] == as.character(level)
    for (q in 1:length(subs)) {
      if (is.na(subs[q])) {
        subs[q] <- F
      }
    }
    x <- x[subs, , drop = F]
    freq <- apply(x, 2, sum)
    subs <- freq > 0
    x <- x[, subs, drop = F]
  }
  x <- as.matrix(x)
  if (dim(x)[1] == 0) {
    result <- array(NA, dim = c(1, 1))
    colnames(result) <- index
    rownames(result) <- "none"
    return(result)
  }
  if (method == "jackknife") {
    thetadiv <- function(x, xdata, index) {
      xdata2 <- xdata[x, 1:ncol(xdata)]
      diversityresult0(xdata2, index = index, method = "pooled")
    }
    if (nrow(x) > 1) {
      result2 <- bootstrap::jackknife(1:nrow(x), thetadiv, 
                                      x, index = index)
      result2$jack.estimate <- mean(as.numeric(result2$jack.values), 
                                    na.rm = T)
    }
    else {
      result2 <- list(jack.values = NA, jack.estimate = NA)
    }
  }
  if (method != "jackknife") {
    result <- diversityresult0(x, index = index, method = method)
  }
  if (method == "mean") {
    result2 <- result
    result <- result2[1]
    result[1] <- mean(result2)
  }
  if (method == "max") {
    result2 <- result
    result <- result2[1]
    result[1] <- max(result2)
  }
  if (method == "sd") {
    result2 <- result
    result <- result2[1]
    result[1] <- sd(result2)
  }
  if (sortit == T && method != "jackknife" && method != "pooled") {
    result <- sort(result)
  }
  if (method != "jackknife") {
    result2 <- round(result, digits = digits)
    result2 <- data.frame(result2)
    colnames(result2) <- index
  }
  if (method == "pooled") {
    rownames(result2) <- "pooled"
  }
  if (method == "mean") {
    rownames(result2) <- "mean"
  }
  if (method == "max") {
    rownames(result2) <- "max"
  }
  if (method == "sd") {
    rownames(result2) <- "sd"
  }
  if (method != "pooled" && method != "jackknife" && method != 
      "mean" && method != "max" && method != "sd") {
    rownames(result2) <- names(result)
  }
  return(result2)
}