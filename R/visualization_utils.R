#' Creates a feature-based data frame to be used by featurePlot or violinPlot
#'
#' This function creates a feature-based data frame to be used by featurePlot
#' or violinPlot
#'
#' @inheritParams documentFun
#' @param feature A gene name or metadata column name.
#'
#' @return A data frame with three columns in which the first two
#' correspond to a dimensionality reduction and the third is the feature score.
#'
#' @keywords internal
#'
pickFeature <- function(scObj, feature, dimred, dims){

    if (length(dims) != 2)
        stop('`dims` must be a vector of size 2.')

    dimred <- dimredName(scObj, dimred)
    df <- as.data.frame(scDimredMat(scObj, dimred)[, dims])

    if (feature %in% rownames(scObj))
        df[, 3] <- scExpMat(scObj, genes=feature)[, 1]

    if (feature %in% metadataNames(scObj)){
        col <- scCol(scObj, feature)
        if (is(col)[1] != 'numeric')
            stop(feature, ' is not a numeric column.')
        df[, 3] <- col
    }

    if(ncol(df) == 2)
        stop('Feature not found: ', feature, '.')

    colnames(df)[3] <- feature
    return(df)
}
