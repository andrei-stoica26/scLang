#' @rdname scDimredMat
#' @export
#'
scDimredMat.default <- function(scObj, dimred)
    stop('Unrecognized input type: `scObj` must be a Seurat or ',
         'SingleCellExpression object.')

#' @rdname scDimredMat
#' @export
#'
scDimredMat.Seurat <- function(scObj, dimred)
{
    reductions <- dimredNames(scObj)
    dimred <- dimredName(scObj, dimred)
    return(as.matrix(Embeddings(scObj, reduction=dimred)))
}

#' @rdname scDimredMat
#' @export
#'
scDimredMat.SingleCellExperiment <- function(scObj, dimred){
    reductions <- dimredNames(scObj)
    dimred <- dimredName(scObj, dimred)
    return(reducedDim(scObj, dimred))
}

#' Extracts the PCA matrix from object.
#'
#' This function extracts the PCA matrix from a Seurat or
#' SingleCellExperiment object.
#'
#' @inheritParams metadataDF
#'
#' @return A PCA matrix.
#'
#' @examples
#' scePath <- system.file('extdata', 'sceObj.qs', package='scLang')
#' sceObj <- qs::qread(scePath)
#' pcaMat <- scPCAMat(sceObj)
#'
#' @export
#'
scPCAMat <- function(scObj)
    return(scDimredMat(scObj, 'pca'))

#' Extracts the UMAP matrix from object.
#'
#' This function extracts the UMAP matrix from a Seurat or
#' SingleCellExperiment object.
#'
#' @inheritParams metadataDF
#'
#' @return A UMAP matrix.
#'
#' @examples
#' scePath <- system.file('extdata', 'sceObj.qs', package='scLang')
#' sceObj <- qs::qread(scePath)
#' umapMat <- scUMAPMat(sceObj)
#'
#' @export
#'
scUMAPMat <- function(scObj)
    return(scDimredMat(scObj, 'umap'))
