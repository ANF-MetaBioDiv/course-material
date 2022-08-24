#' Transfer row.names values to a new column
#'
#'
#' @return Input data.frame with an extra column contains
#' @param rich A dataframe with alpha diversity indices as columns
#' and samples as rows as the output of pyloseq::estimate_richness()
#' @param file Output file name

df_export <- function(input_table, new_rn = "variable") {
    tmp <- data.frame(row.names(input_table),
        input_table,
        row.names = NULL)
    names(tmp)[1] <- new_rn
    return(tmp)
}
