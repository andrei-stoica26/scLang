#' Extract feature data to be used by \code{featurePlot} or \code{violinPlot}
#'
#' This function extracts feature data to be used by \code{featurePlot}
#' or \code{violinPlot}.
#'
#' @inheritParams documentFun
#' @param feature A gene name or metadata column name.
#'
#' @return A numeric vector.
#'
#' @keywords internal
#'
pickFeature <- function(scObj, feature){
    if (feature %in% rownames(scObj))
        return(scExpMat(scObj, genes=feature)[, 1])

    if (feature %in% metadataNames(scObj)){
        col <- scCol(scObj, feature)
        if (!is.numeric(col))
            stop(feature, ' is not a numeric column.')
        return(col)
    }

    stop('Feature not found: ', feature, '.')
}
