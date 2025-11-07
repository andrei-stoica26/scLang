#' @importFrom dplyr count
#'
NULL

#' Extract count information from single-cell expression object column
#'
#' This function extracts count information from the column of a Seurat or
#' SingleCellExperiment object.
#'
#' @inheritParams metadataDF
#' @param col Column as string.
#'
#' @return A frequency vector with the unique column values as names.
#'
#' @examples
#' scePath <- system.file('extdata', 'sceObj.qs', package='scLang')
#' sceObj <- qs::qread(scePath)
#' scColCounts(sceObj, 'Mutation_Status')
#'
#' @export
#'
scColCounts <- function(scObj, col='orig.ident'){
    df <- dplyr::count(metadataDF(scObj), .data[[col]])
    v <- setNames(df[, 2], as.factor(df[, 1]))
    return(v)
}

#' Extract count information from Seurat column
#'
#' This function extracts count information from Seurat column.
#'
#' @inheritParams metadataDF
#' @param col1 Column as string.
#' @param col2 Column as string.
#'
#' @return A data frame listing the counts of all combinations of pairs from
#' two categorical columns.
#'
#' @examples
#' scePath <- system.file('extdata', 'sceObj.qs', package='scLang')
#' sceObj <- qs::qread(scePath)
#' scColPairCounts(sceObj, 'Mutation_Status', 'Cell_Cycle')
#'
#' @export
#'
scColPairCounts <- function(scObj, col1='seurat_clusters', col2='orig.ident')
    return(dplyr::count(metadataDF(scObj), .data[[col1]], .data[[col2]]))
