#' @importFrom dplyr count
#'
NULL

#' Extract per-group counts from the column of a single-cell expression object
#'
#' This function extracts per-group counts from the column of single-cell
#' expression object.
#'
#' @inheritParams metadataDF
#' @param col Column as string.
#'
#' @return A frequency vector with the unique column values as names.
#'
#' @examples
#' scePath <- system.file('extdata', 'sceObj.qs2', package='scLang')
#' sceObj <- qs2::qs_read(scePath)
#' scColCounts(sceObj, 'Mutation_Status')
#'
#' @export
#'
scColCounts <- function(scObj, col='orig.ident'){
    df <- dplyr::count(metadataDF(scObj), .data[[col]])
    v <- setNames(df[, 2], as.factor(df[, 1]))
    return(v)
}

#' Extract counts from two columns of a single-cell expression
#' object
#'
#' This function extracts count information from two columns of a single-cell
#' expression object.
#'
#' @inheritParams metadataDF
#' @param col1 Column as string.
#' @param col2 Column as string.
#'
#' @return A data frame listing the counts of all combinations of pairs from
#' two categorical columns.
#'
#' @examples
#' scePath <- system.file('extdata', 'sceObj.qs2', package='scLang')
#' sceObj <- qs2::qs_read(scePath)
#' scColPairCounts(sceObj, 'Mutation_Status', 'Cell_Cycle')
#'
#' @export
#'
scColPairCounts <- function(scObj, col1='seurat_clusters', col2='orig.ident')
    return(dplyr::count(metadataDF(scObj), .data[[col1]], .data[[col2]]))

#' Extract percentages from two columns of a single-cell expression
#' object
#'
#' This function extracts percentage information from two columns of a single-cell
#' expression object.
#'
#' @inheritParams scColPairCounts
#'
#' @return A data frame listing the percentages of all combinations of pairs
#' from two categorical columns.
#'
#' @examples
#' scePath <- system.file('extdata', 'sceObj.qs2', package='scLang')
#' sceObj <- qs2::qs_read(scePath)
#' scColPairPercs(sceObj, 'Mutation_Status', 'Cell_Cycle')
#'
#' @export
#'
scColPairPercs <- function(scObj, col1, col2, sigDigits = 2){
    df <- scColPairCounts(scObj, col1, col2)
    nRep <- length(unique(scCol(scObj, col2)))
    totals <- as.numeric(unlist(lapply(scColCounts(scObj, col1),
                                       function(x) rep(x, nRep))))
    df$perc <- round(df$n / totals * 100, sigDigits)
    return(df)
}
