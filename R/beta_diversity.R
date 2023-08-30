
run_ancombc_patched <- function(
    ps,
    group,
    confounders = character(0),
    contrast = NULL,
    taxa_rank = "all",
    transform = c("identity", "log10", "log10p"),
    norm = "none",
    norm_para = list(),
    p_adjust = c(
        "none", "fdr", "bonferroni", "holm",
        "hochberg", "hommel", "BH", "BY"
    ),
    prv_cut = 0.1,
    lib_cut = 0,
    struc_zero = FALSE,
    neg_lb = FALSE,
    tol = 1e-05,
    max_iter = 100,
    conserve = FALSE,
    pvalue_cutoff = 0.05) {

    stopifnot(inherits(ps, "phyloseq"))

    ps <- microbiomeMarker:::check_rank_names(ps) |>
        microbiomeMarker:::check_taxa_rank(taxa_rank)

    if (length(confounders)) {
        confounders <- microbiomeMarker:::check_confounder(ps, group, confounders)
    }

    # if it contains missing values for any
    # variable specified in the formula, the corresponding sampling fraction
    # estimate for this sample will return NA since the sampling fraction is
    # not estimable with the presence of missing values.
    # remove this samples
    fml_char <- ifelse(
        length(confounders),
        paste(c(confounders, group), collapse = " + "),
        group
    )

    # fml_char <- paste(c(confounders, group), collapse = " + ")
    # fml <- stats::as.formula(paste("~", fml_char))
    # vars_fml <- all.vars(fml)

    for (var in c(confounders, group)) {
        ps <- microbiomeMarker:::remove_na_samples(ps, var)
    }

    # check whether group is valid, write a function
    meta <- phyloseq::sample_data(ps)

    groups <- meta[[group]]
    groups <- make.names(groups)
    if (!is.null(contrast)) {
        contrast <- make.names(contrast)
    }
    if (!is.factor(groups)) {
        groups <- factor(groups)
    }
    groups <- microbiomeMarker:::set_lvl(groups, contrast)
    phyloseq::sample_data(ps)[[group]] <- groups
    lvl <- levels(groups)
    n_lvl <- length(lvl)

    contrast <- microbiomeMarker:::check_contrast(contrast)

    transform <- match.arg(transform, c("identity", "log10", "log10p"))
    p_adjust <- match.arg(
        p_adjust,
        c(
            "none", "fdr", "bonferroni", "holm",
            "hochberg", "hommel", "BH", "BY"
        )
    )

    # set the reference level for pair-wise comparison from mutliple groups
    # if (!is.null(contrast) && n_lvl > 2) {
    #     groups <- relevel(groups, ref = contrast[1])
    #     sample_data(ps)[[group]] <- groups
    # }


    # preprocess phyloseq object
    ps <- microbiomeMarker:::preprocess_ps(ps)
    ps <- microbiomeMarker::transform_abundances(ps, transform = transform)

    # normalize the data
    norm_para <- c(norm_para, method = norm, object = list(ps))
    ps_normed <- do.call(microbiomeMarker::normalize, norm_para)
    ps_summarized <- microbiomeMarker:::pre_ps_taxa_rank(ps_normed, taxa_rank)

    global <- ifelse(n_lvl > 2, TRUE, FALSE)
    # ancombc differential abundance analysis

    if (taxa_rank == "all") {
        ancombc_taxa_rank <- phyloseq::rank_names(ps_summarized)[1]
    } else if(taxa_rank == "none") {
        ancombc_taxa_rank <- NULL
    } else {
        ancombc_taxa_rank <- taxa_rank
    }

    ancombc_out <- ANCOMBC::ancombc(
        ps_summarized,
        tax_level = ancombc_taxa_rank,
        formula = fml_char,
        p_adj_method = p_adjust,
        prv_cut = prv_cut,
        lib_cut = lib_cut,
        group = group,
        struc_zero = struc_zero,
        neg_lb = neg_lb,
        tol = tol,
        max_iter = max_iter,
        conserve = conserve,
        alpha = pvalue_cutoff,
        global = global
    )

    # multiple-group comparison will be performed while the group
    # variable has > 2 levels
    keep_var <- c("W", "p_val", "q_val", "diff_abn")
    if (n_lvl > 2) {
        # ANCOM-BC global test to determine taxa that are differentially
        # abundant between three or more groups of multiple samples.
        # global result to marker_table
        if (is.null(contrast)) {
            mtab <- ancombc_out$res_global
        } else {
            exp_lvl <- paste0(group, contrast[2])
            ancombc_res <- ancombc_out$res
            mtab <- lapply(keep_var, function(x) ancombc_res[[x]][exp_lvl])
            mtab <- do.call(cbind, mtab)
        }
    } else {
        ancombc_out_res <- ancombc_out$res
        # drop intercept
        ancombc_out_res <- lapply(
            ancombc_out_res,
            function(x) data.frame(x, row.names = "taxon")[-1])
        mtab <- do.call(
            cbind,
            ancombc_out_res[c("W", "p_val", "q_val", "diff_abn")]
        )
    }
    names(mtab) <- keep_var

    # determine enrich group based on coefficients
    # drop intercept
    cf <- data.frame(ancombc_out$res$lfc, row.names = "taxon")[-1]
    if (n_lvl > 2) {
        if (!is.null(contrast)) {
            cf <- cf[exp_lvl]
            enrich_group <- ifelse(cf[[1]] > 0, contrast[2], contrast[1])
        } else {
            cf <- cbind(0, cf)
            enrich_group <- lvl[apply(cf, 1, which.max)]
        }
    } else {
        enrich_group <- ifelse(cf[[1]] > 0, lvl[2], lvl[1])
    }

    # # enriched group
    enrich_abd <- microbiomeMarker:::get_ancombc_enrich_group(ps_summarized, ancombc_out, group)
    norm_abd <- enrich_abd$abd

    mtab <- cbind(feature = row.names(mtab), mtab, enrich_group)
    mtab_sig <- mtab[mtab$diff_abn, ]
    mtab_sig <- mtab_sig[c("feature", "enrich_group", "W", "p_val", "q_val")]
    names(mtab_sig) <- c("feature", "enrich_group", "ef_W", "pvalue", "padj")
    marker <- microbiomeMarker:::return_marker(mtab_sig, mtab)
    marker <- microbiomeMarker::microbiomeMarker(
        marker_table = marker,
        norm_method = microbiomeMarker:::get_norm_method(norm),
        diff_method = "ancombc",
        sam_data = phyloseq::sample_data(ps_summarized),
        otu_table = phyloseq::otu_table(norm_abd, taxa_are_rows = TRUE),
        tax_table = phyloseq::tax_table(ps_summarized)
    )

    marker
}