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
#' scColCounts(sceObj, 'Cluster')
#'
#' @export
#'
scColCounts <- function(scObj, col){
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
#' scColPairCounts(sceObj, 'Cluster', 'Donor')
#'
#' @export
#'
scColPairCounts <- function(scObj, col1, col2)
    return(dplyr::count(metadataDF(scObj), .data[[col1]], .data[[col2]]))

#' Extract percentages from two columns of a single-cell expression
#' object
#'
#' This function extracts percentage information from two columns of a
#' single-cell expression object. For each i x j combination with i taken
#' from column 1 and j taken from column 2, the function reports the percentage
#' that j contributes to all combinations involving i.
#'
#' @inheritParams scColPairCounts
#' @param sigDigits Number of significant digits.
#'
#' @return A data frame listing the percentages of all combinations of pairs
#' from two categorical columns.
#'
#' @examples
#' scePath <- system.file('extdata', 'sceObj.qs2', package='scLang')
#' sceObj <- qs2::qs_read(scePath)
#' scColPairPercs(sceObj, 'Cluster', 'Donor')
#'
#' @export
#'
scColPairPercs <- function(scObj, col1, col2, sigDigits = 2){
    df <- scColPairCounts(scObj, col1, col2)
    totals <- vapply(unique(df[, 1]), function(x)
        sum(df[df[, 1] == x, ]$n), integer(1))
    df$perc <- round(df[, 3] / totals[df[, 1]] * 100, sigDigits)
    return(df)
}

#' Extract observed-over-expected ratios from two columns of a single-cell
#' expression object
#'
#' This function extracts observed-over-expected ratios from two columns of a
#' single-cell expression object.
#'
#' @inheritParams scColPairCounts
#' @param sigDigits Number of significant digits.
#'
#' @return A data frame listing the observed-over-expected ratios of all
#' combinations of pairs from two categorical columns.
#'
#' @examples
#' scePath <- system.file('extdata', 'sceObj.qs2', package='scLang')
#' sceObj <- qs2::qs_read(scePath)
#' scColPairRatios(sceObj, 'Cluster', 'Donor')
#'
#' @export
#'
scColPairRatios <- function (scObj, col1, col2, sigDigits = 2){
    df <- scColPairPercs(scObj, col1, col2)
    totals <- vapply(unique(df[, 2]), function(x)
        sum(df[df[, 2] == x, ]$n), integer(1))
    total <- sum(totals)
    df$totalPerc <- round(totals[df[, 2]] / total * 100, sigDigits)
    df$ratio <- df$perc / df$totalPerc
    return(df)
}
