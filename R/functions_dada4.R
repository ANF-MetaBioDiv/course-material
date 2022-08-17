################### MES FONCTIONS POUR DADA2 ###################################

####################Qualité profile avec plotQualityProfile ##################
#retourne un fichier pdf 

qualityprofile <- function(fileFw,fileRev,outfile)
{
  pdf(outfile, onefile = TRUE)
  for (i in 1:length(fileFw))
  {
    f <- dada2::plotQualityProfile(fileFw[[i]])
    r <- dada2::plotQualityProfile(fileRev[[i]])
    print(f)
    print(r)
  }
  dev.off()
}
###############################################################################

############# COMPTAGE des amorces FWD et REV dans les Fichiers seq ############
######### 1-Fonction pour mettre primers en Fwd, rev, et Revcomplement : toutes les orientations possibles
allOrients <- function(primer)
{
  #require(Biostrings)
  dna <- Biostrings::DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Reverse = reverse(dna), RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
###############################################################################
########## Fonction Compter hit de primers #########################
primerHits <- function(primer, fn) 
  {
  #require(ShortRead)
    # Counts number of reads in which the primer is found
  nhits <- Biostrings::vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
  }

######################################################################
############ Fonction: boucle de comptage sur tous les fichiers ###### 
count<- function(infile,whichprimer,out) 
{
  for (i in 1:length(infile))
  {
    nom <-basename(infile[[i]]) #recuperer qu'une partie du chemin  pour en faire le nom
    print(paste("Recherche Amorces dans:",nom))
    #boucle pour comptage du primer  & ecriture en continue (append) dans Fichier out
    all <- sapply(whichprimer, primerHits, fn = infile[[i]])
    write.table(c(nom,all),out, append=TRUE,col.names = F)
  }
}

#######################################################################

#######################Fonction Fastqc Plot ##########################################
##necessite package fastqcr
#Repertoire contenant vos fichiers bruts: fq.dir et qc.dir : repertoire pour les resultats

fastqcanalyse <- function (path,path1)
{
  fastqc(fq.dir = path, qc.dir = path1, threads = 4)
  #mettre dans variable qc.files la liste des chemins ou sont les fichiers zip issues du resultat fastqc
  qc.dir=path1
  qc.files <- list.files(path1, full.names = TRUE, pattern = ".zip")
  #recuperer le nom sample 
  sample.names_fastqc <- sapply(strsplit(basename(qc.files), "_"), `[`, 1)
  #lire la collection de ces fichiers et la mettre ds variable qc
  qc <- qc_read_collection(qc.files, sample_names=sample.names_fastqc)
  
  #creer un fichier html report fastqc 
  qc <- qc_aggregate(qc.dir, progressbar = TRUE)
  qc_report(qc.dir, result.file = path2, experiment = "mydata")
  
  #on peut maintenant executer les plots pour tous les fichiers et diriger vers un pdf
  #pdf("mygraph.pdf")     # création de fichier pdf 
  #qc_plot_collection(qc, "Per base sequence quality")
  #qc_plot_collection(qc, "Per sequence quality scores")
  #qc_plot_collection(qc, "Sequence Length Distribution")
  #dev.off() # pour arrêter la redirection du graphique 
  
}
###########Verifier Normalité des donnees de tous les Indices alpha pour les tests stats
### Prend en entree le fichier de resultat de la fonction phyloseq "estimate_richness"

indices_normality <- function (rich)
{
  ### p-value < 0.05 means data failed normality test
  columns <- names(rich)
  pdf("./Distribution_Indices_normality.pdf", onefile = FALSE)
  
  par(mfrow=c(3,3))
  
  for (i in columns) 
  {
    shap <- shapiro.test(rich[,i])
    qqnorm(rich[,i],main=i, sub=shap$p.value);qqline(rich[,i])
   # hist(rich[,i], main=i, sub=shap$p.value, xlab="", breaks=10)
    #lines(density(rich[,i]), col = "red",lwd=2)
  }
  dev.off()
}


########## Fonction recuperation des sequences des ASVs a a partir d'une liste
################################################
getseq <- function(listID, objphylo,file)
{
  seqq<- (DNAStringSet())
  for (i in listID)
  {
    nom<- taxa_names(objphylo@refseq[i]);
    seqq[i]<- (DNAStringSet(objphylo@refseq[i]));
    writeXStringSet(seqq,file,format="fasta")
  }
 
}
##################################################
#' @title Count number of OTUs by taxonomic rank for each sample.
#' @details This function will split phyloseq object by the specified taxonomic rank
#' @param physeq A phyloseq-class object
#' @param TaxRank Name of the taxonomic rank
#' @param relative Logical, return relative number of OTUs
#' @param add_meta_data Logical, add sample metadata to the resulting table
#'
#' @return Data frame with OTU counts (columns: Sample, TaxRank, N.OTU + optional sample metadata).
#' @export
#'
#' @examples
#' data(GlobalPatterns)
#' phylum_counts <- phyloseq_ntaxa_by_tax(GlobalPatterns, TaxRank = "Phylum")
#' head(phylum_counts)
#'
phyloseq_ntaxa_by_tax <- function(physeq, TaxRank = "Phylum", relative = F, add_meta_data = T){
  
  ## Melt phyloseq data object into large data.frame
  mm <- phyloseq::psmelt(physeq)
  
  ## Count number of OTUs for each sample
  count_otus <- function(z, TaxRank = TaxRank, relative = relative){
    
    ## Remove zero-OTUs
    zero_otu <- z$Abundance > 0
    if(any(zero_otu)){ z <- z[zero_otu, ] }
    
    ## Count number of OTUs per taxonomic rank selected
    rez <- as.data.frame(table(z[, TaxRank]), stringsAsFactors = F)
    colnames(rez) <- c(TaxRank, "N.OTU")
    
    ## Transform to relative abundance
    if(relative == TRUE){
      rez$N.OTU <- with(rez, N.OTU / sum(N.OTU) )
    }
    
    return(rez)
  }
  
  ## Count number of OTUs for each sample
  res <- plyr::ddply(.data = mm, .variables = "Sample", .fun = count_otus, TaxRank = TaxRank, relative = relative)
  
  ## Add meta-data
  if(add_meta_data == TRUE){
    if(!is.null(phyloseq::sample_data(physeq, errorIfNULL = FALSE))){  # add only if metadata is present
      ## Extract meta-data
      metad <- data.frame(Sample = phyloseq::sample_names(physeq), phyloseq::sample_data(physeq))
      
      ## Extract column names
      # main_cols <- c("OTU", "Sample", "Abundance", rank_names(x))          # 'standard' columns
      # meta_cols <- colnames(mm)[ which(!colnames(mm) %in% main_cols) ]     # meta-data columns
      main_cols <- c("Sample")                                               # 'standard' columns
      meta_cols <- colnames(metad)[ which(!colnames(metad) %in% main_cols) ] # meta-data columns
      
      res <- cbind(res, metad[match(x = res$Sample, table = metad$Sample), meta_cols] )
    }
  }
  
  rownames(res) <- NULL
  return(res)
}
################################################
# @ Thomas W. Battaglia

#' Estimate phylogenetic diversity from phyloseq object pd_whole_tree
#'
#' Estimate the Faiths phylogenetic diverstiy from an OTU table and phylogenetic tree
#'
#' @param phylo A phyloseq object with an OTU table and phylogenetic tree slot.
#' @return A data.frame of phylogenetic diversity metrics.
#' @export
# Estimate PD-whole tree
estimate_pd <- function(phylo){
  
  # Error if input is not of class phylo
  if(class(phylo) != "phyloseq"){
    stop("Input file is not of class 'phyloseq'.")
  }
  
  # Error if no class phy_tree
  if(!(.hasSlot(phylo, "phy_tree"))){
    stop("Could not find tree slot in phylo object.")
  }
  
  # Transpose if needed
  # Adapted from phyloseq/vegan import
  OTU <- phyloseq::otu_table(phylo)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  
  # Get matrix version of OTU table
  otutable <- as(OTU, "matrix")
  
  # Get phylogenetic tree from pyloseq object
  tree <- phyloseq::phy_tree(phylo)
  
  # Print status message
  message("Calculating Faiths PD-index...")
  
  # If object is greater than 10mb, then print status message
  if(object.size(otutable) > 10000000){
    message("This is a large object, it may take awhile...")
  }
  
  # Calculate Faith's PD-index
  pdtable <- picante::pd(otutable, tree, include.root = F)
  
  # Return data frame of results
  return(pdtable)
}
#############################################################@@

ggrare <- function(physeq, step = 10, label = NULL, color = NULL, plot = TRUE, parallel = FALSE, se = TRUE) {
  ## Args:
  ## - physeq: phyloseq class object, from which abundance data are extracted
  ## - step: Step size for sample size in rarefaction curves
  ## - label: Default `NULL`. Character string. The name of the variable
  ##          to map to text labels on the plot. Similar to color option
  ##          but for plotting text.
  ## - color: (Optional). Default ‘NULL’. Character string. The name of the
  ##          variable to map to colors in the plot. This can be a sample
  ##          variable (among the set returned by
  ##          ‘sample_variables(physeq)’ ) or taxonomic rank (among the set
  ##          returned by ‘rank_names(physeq)’).
  ##
  ##          Finally, The color scheme is chosen automatically by
  ##          ‘link{ggplot}’, but it can be modified afterward with an
  ##          additional layer using ‘scale_color_manual’.
  ## - color: Default `NULL`. Character string. The name of the variable
  ##          to map to text labels on the plot. Similar to color option
  ##          but for plotting text.
  ## - plot:  Logical, should the graphic be plotted.
  ## - parallel: should rarefaction be parallelized (using parallel framework)
  ## - se:    Default TRUE. Logical. Should standard errors be computed. 
  ## require vegan
  x <- as(otu_table(physeq), "matrix")
  if (taxa_are_rows(physeq)) { x <- t(x) }
  
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
    y <- rarefy(x[i, ,drop = FALSE], n, se = se)
    if (nrow(y) != 1) {
      rownames(y) <- c(".S", ".se")
      return(data.frame(t(y), Size = n, Sample = rownames(x)[i]))
    } else {
      return(data.frame(.S = y[1, ], Size = n, Sample = rownames(x)[i]))
    }
  }
  if (parallel) {
    out <- mclapply(seq_len(nr), rarefun, mc.preschedule = FALSE)
  } else {
    out <- lapply(seq_len(nr), rarefun)
  }
  df <- do.call(rbind, out)
  
  ## Get sample data 
  if (!is.null(sample_data(physeq, FALSE))) {
    sdf <- as(sample_data(physeq), "data.frame")
    sdf$Sample <- rownames(sdf)
    data <- merge(df, sdf, by = "Sample")
    labels <- data.frame(x = tot, y = S, Sample = rownames(x))
    labels <- merge(labels, sdf, by = "Sample")
  }
  
  ## Add, any custom-supplied plot-mapped variables
  if( length(color) > 1 ){
    data$color <- color
    names(data)[names(data)=="color"] <- deparse(substitute(color))
    color <- deparse(substitute(color))
  }
  if( length(label) > 1 ){
    labels$label <- label
    names(labels)[names(labels)=="label"] <- deparse(substitute(label))
    label <- deparse(substitute(label))
  }
  
  p <- ggplot(data = data, aes_string(x = "Size", y = ".S", group = "Sample", color = color))
  p <- p + labs(x = "Sample Size", y = "Species Richness")
  if (!is.null(label)) {
    p <- p + geom_text(data = labels, aes_string(x = "x", y = "y", label = label, color = color),
                       size = 4, hjust = 0)
  }
  p <- p + geom_line()
  if (se) { ## add standard error if available
    p <- p + geom_ribbon(aes_string(ymin = ".S - .se", ymax = ".S + .se", color = NULL, fill = color), alpha = 0.2)
  }
  if (plot) {
    plot(p)
  }
  invisible(p)
}

############################

#@name ps_extra-accessors
#' @title Extract elements from ps_extra class
#'
#' @description
#' - `ps_get`     returns phyloseq
#' - `info_get`   returns ps_extra_info object
#' - `dist_get`   returns distance matrix (or NULL)
#' - `ord_get`    returns ordination object (or NULL)
#' - `perm_get`   returns adonis2() permanova model (or NULL)
#' - `bdisp_get`  returns results of betadisper() (or NULL)
#' - `otu_get`    returns phyloseq otu_table matrix with taxa as columns
#' - `tt_get`     returns phyloseq tax_table
#' - `samdat_tbl` returns phyloseq sample_data as a tibble,
#' with sample_names as new first column called .sample_name
#'
#' @param ps_extra ps_extra class object
#'
#' @return element of ps_extra class object (or NULL)
#' @export
#'
#' @examples
#' data("dietswap", package = "microbiome")
#' psx <- tax_transform(dietswap, "identity", rank = "Genus")
#' psx
#'
#' ps_get(psx)
#' info_get(psx)
#'
#' dist_get(psx) # this ps_extra has no dist_calc result
#' ord_get(psx) # this ps_extra has no ord_calc result
#' perm_get(psx) # this ps_extra has no dist_permanova result
#' bdisp_get(psx) # this ps_extra has no dist_bdisp result
#'
#' # these can be returned from phyloseq objects too
#' otu_get(psx)[1:6, 1:4]
#' tt_get(psx) %>% head()
#' samdat_tbl(psx) %>% head()
#' @export
#' @rdname ps_extra-accessors
ps_get <- function(ps_extra) {
  if (inherits(ps_extra, "ps_extra")) {
    return(ps_extra[["ps"]])
  }
  if (inherits(ps_extra, "phyloseq")) {
    return(ps_extra)
  }
  stop(
    'class of argument should be "ps_extra" or "phyloseq"',
    "\nargument is class: ", paste(class(ps_extra), collapse = " ")
  )
}
#' @rdname ps_extra-accessors
#' @export
dist_get <- function(ps_extra) {
  stopifnot(inherits(ps_extra, "ps_extra"))
  ps_extra[["dist"]]
}
#' @rdname ps_extra-accessors
#' @export
ord_get <- function(ps_extra) {
  stopifnot(inherits(ps_extra, "ps_extra"))
  ps_extra[["ord"]]
}
#' @rdname ps_extra-accessors
#' @export
info_get <- function(ps_extra) {
  if (inherits(ps_extra, "ps_extra")) {
    return(ps_extra[["info"]])
  } else if (methods::is(ps_extra, "phyloseq")) {
    return(new_ps_extra_info())
  } else {
    stop(
      "info_get can only get info from a 'ps_extra' class object.",
      "\nargument is class: ", paste(class(ps_extra), collapse = " ")
    )
  }
}
#' @rdname ps_extra-accessors
#' @export
perm_get <- function(ps_extra) {
  stopifnot(inherits(ps_extra, "ps_extra"))
  ps_extra[["permanova"]]
}
#' @rdname ps_extra-accessors
#' @export
bdisp_get <- function(ps_extra) {
  stopifnot(inherits(ps_extra, "ps_extra"))
  ps_extra[["bdisp"]]
}


#' @param data phyloseq or ps_extra
# @return phyloseq otu_table matrix with taxa as columns
#'
#' @param taxa subset of taxa to return, NA for all (default)
#' @param samples subset of samples to return, NA for all (default)
#' @param counts should otu_get ensure it returns counts? if present in object
#'
#' @rdname ps_extra-accessors
#' @export
otu_get <- function(data, taxa = NA, samples = NA, counts = FALSE) {
  # get otu_table from object
  if (methods::is(data, "otu_table")) {
    if (isTRUE(counts)) warning("data is otu_table: ignoring `counts = TRUE`")
    otu <- data
  } else {
    if (isTRUE(counts)) ps <- ps_counts(data)
    if (!isTRUE(counts)) ps <- ps_get(data)
    otu <- phyloseq::otu_table(ps)
  }
  if (phyloseq::taxa_are_rows(otu)) otu <- phyloseq::t(otu)
  
  # subset samples and or taxa if requested
  if (!identical(taxa, NA) || !identical(samples, NA)) {
    otu <- otu[samples, taxa, drop = FALSE]
  }
  return(otu)
}

#' @rdname ps_extra-accessors
#' @export
tt_get <- function(data) {
  if (!methods::is(data, "taxonomyTable")) {
    ps <- ps_get(data)
    tt <- phyloseq::tax_table(ps)
  } else {
    tt <- data
  }
  return(tt)
}

#' @param data phyloseq or ps_extra
# @return phyloseq sample_data as a tibble,
# with sample_names as new first column called .sample_name
#' @param sample_names_col
#' name of column where sample_names are put.
#' if NA, return data.frame with rownames (sample_names)
#' @rdname ps_extra-accessors
#' @export
samdat_tbl <- function(data, sample_names_col = ".sample_name") {
  if (inherits(data, "ps_extra")) data <- ps_get(data)
  if (methods::is(data, "phyloseq") || methods::is(data, "sample_data")) {
    df <- samdatAsDataframe(data)
  } else {
    stop(
      "data must be of class 'phyloseq', 'ps_extra', or 'sample_data', not: ",
      paste(class(data), collapse = " ")
    )
  }
  if (identical(sample_names_col, NA)) {
    return(df)
  } else {
    df <- tibble::rownames_to_column(df, var = sample_names_col)
    return(tibble::as_tibble(df))
  }
}

# internal helper that get phyloseq sample_data as plain dataframe
# without changing invalid colnames (like microbiome::meta does)
# or losing rownames / sample_names (like data.frame() with defaults does)
samdatAsDataframe <- function(ps) {
  samdat <- phyloseq::sample_data(ps)
  df <- data.frame(samdat, check.names = FALSE)
  return(df)
}

# get phyloseq with counts if available
ps_counts <- function(data, warn = TRUE) {
  # always get ps, regardless of ps_extra or phyloseq data or counts presence
  ps <- ps_get(data)
  # checking names of a ps will return NULL (and x %in% NULL returns FALSE)
  if ("counts" %in% names(data)) {
    # get counts and use them if they exist,
    # and check regardless if otutab returned will be counts
    counts <- data[["counts"]]
    # maintain existing taxa_are_rows status for consistency
    if (phyloseq::taxa_are_rows(ps)) counts <- phyloseq::t(counts)
    phyloseq::otu_table(ps) <- counts
  }
  if (!isFALSE(warn)) {
    mess <- paste0(
      "otu_table of counts is NOT available!\n",
      "Available otu_table contains non-zero values that are less than 1"
    )
    # now check ps otu_table is counts
    test_matrix <- unclass(otu_get(ps))
    if (any(test_matrix < 1 & test_matrix != 0)) {
      if (isTRUE(warn)) {
        warning(mess)
      } else if (identical(warn, "error")) {
        stop(mess)
      } else {
        stop("warn argument value is invalid: should be T, F or 'error'")
      }
    }
  }
  return(ps)
}

# ps_extra methods for phyloseq accessors -------------------------------------
methods::setOldClass("ps_extra")

# methods::setGeneric("otu_table", def = phyloseq::otu_table)
methods::setMethod(
  f = phyloseq::otu_table, signature = c(object = "ps_extra"),
  definition = function(object) {
    warning(
      "Using otu_table() with ps_extra objects is not recommended.\n\t",
      "Use otu_get() instead, which always returns taxa as columns."
    )
    return(phyloseq::otu_table(ps_get(object)))
  }
)

methods::setMethod(
  f = phyloseq::sample_data, signature = c(object = "ps_extra"),
  definition = function(object) phyloseq::sample_data(ps_get(object))
)

methods::setMethod(
  f = phyloseq::tax_table, signature = c(object = "ps_extra"),
  definition = function(object) {
    return(phyloseq::tax_table(ps_get(object)))
  }
)

methods::setMethod(
  f = phyloseq::sample_names, signature = c(physeq = "ps_extra"),
  definition = function(physeq) phyloseq::sample_names(ps_get(physeq))
)

methods::setMethod(
  f = phyloseq::taxa_names, signature = c(physeq = "ps_extra"),
  definition = function(physeq) phyloseq::taxa_names(ps_get(physeq))
)

methods::setMethod(
  f = phyloseq::phy_tree, signature = c(physeq = "ps_extra"),
  definition = function(physeq) phyloseq::phy_tree(ps_get(physeq))
)

methods::setMethod(
  f = phyloseq::refseq, signature = c(physeq = "ps_extra"),
  definition = function(physeq) phyloseq::refseq(ps_get(physeq))
)

# rank names is not a generic in phyloseq
methods::setGeneric(name = "rank_names", def = phyloseq::rank_names)
methods::setMethod(
  f = "rank_names", signature = c(physeq = "ps_extra"),
  definition = function(physeq) phyloseq::rank_names(ps_get(physeq))
)

########## PS_reorder
#' Set order of samples in phyloseq object
#'
#' @description
#' Manually set order of samples by specifying samples names in desired order.
#'
#' @details
#' Ordering of samples in a phyloseq is controlled from the otu_table slot!
#'
#' @param ps phyloseq
#' @param sample_order names or current numerical indices of samples in desired order
#'
#' @return phyloseq
#' @export
#'
#' @seealso \code{\link{ps_arrange}} for arranging samples by sample_data variables (or otu_table)
#' @seealso \code{\link{ps_seriate}} for arranging samples by microbiome similarity
#' @seealso \code{\link{ps_filter}} for keeping only some samples, based on sample_data
#'
#' @examples
#' library(phyloseq)
#' data("dietswap", package = "microbiome")
#'
#' dietswap %>%
#'   sample_data() %>%
#'   head(8)
#'
#' new_order <- rev(sample_names(dietswap))
#' dietswap %>%
#'   ps_reorder(new_order) %>%
#'   sample_data() %>%
#'   head(8)
#'
#' # random ordering with numbers
#' set.seed(1000)
#' random_order <- sample(1:nsamples(dietswap))
#' dietswap %>%
#'   ps_reorder(random_order) %>%
#'   sample_data() %>%
#'   head(8)
ps_reorder <- function(ps, sample_order) {
  ps <- ps_get(ps)
  otu <- phyloseq::otu_table(ps)
  if (phyloseq::taxa_are_rows(ps)) {
    otu <- otu[, sample_order]
  } else {
    otu <- otu[sample_order, ]
  }
  phyloseq::otu_table(ps) <- otu
  
  return(ps)
}

#########Good's coverage 
#Calculate Good's Coverage
#'
#' Calculates Good's coverage from a community data matrix with samples as rows and OTUs as columns.
#'
#' @param com a vegan compatible community data matrix.
#'
#' @return A table with the headings number of singletons, number of sequences, and Good's coverage for each sample in rows.
#' @export
#'
#' @references Good, I. J. 1953. The Population Frequencies of Species and the Estimation of Population Parameters. Biometrika 40:237-264.
#'
#' @examples
#'
goods <-
  function(com){
    no.seqs <- rowSums(com)
    sing <- com==1
    no.sing <- apply(sing, 1, sum)
    goods <- 100*(1-no.sing/no.seqs)
    goods.sum <- cbind(no.sing, no.seqs, goods)
    goods.sum <- as.data.frame(goods.sum)
    return(goods.sum)
  }



