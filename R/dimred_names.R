#' Correct the case of the input dimensionality reduction name if necessary
#'
#' This function corrects, if necessary, the case of the input dimensionality
#' reduction name. The dimensionality reduction must be available in the object
#' whether in upper or in lower case. Otherwise, the function will return
#' an error.
#'
#' @inheritParams scDimredMat
#'
#' @return The dimensionality reduction name.
#'
#' @noRd
#'
dimredName <- function(scObj, dimred){
    lowerDimred <- tolower(dimred)
    upperDimred <- toupper(dimred)
    dimName <- intersect(dimredNames(scObj), c(dimred, lowerDimred,
                                               upperDimred))
    if (!length(dimName))
        stop('No `', lowerDimred, '` or `', upperDimred,
             '` reduction was found in the object.')
    if (length(dimName) > 1){
        warning('Both `', lowerDimred, '` and ', upperDimred, '` reductions were ',
                'found. Returning ', dimred, '.')
        return(dimred)
    }
    return(dimName)
}

###############################################################################

#' @rdname dimredNames
#' @export
#'
dimredNames.default <- function(scObj)
    stop('Unrecognized input type: `scObj` must be a Seurat, ',
         'SingleCellExpression, matrix or dgCMatrix object.')

#' @rdname dimredNames
#' @export
#'
dimredNames.Seurat <- function(scObj)
    return(Reductions(scObj))


#' @rdname dimredNames
#' @export
#'
dimredNames.SingleCellExperiment <- function(scObj)
    return(names(reducedDims(scObj)))
