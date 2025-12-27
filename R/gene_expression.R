#' @importFrom rlang .data
#' @importFrom S4Vectors DataFrame
#' @importFrom SeuratObject Embeddings LayerData Layers Reductions
#' @importFrom SingleCellExperiment reducedDim reducedDims
#' @importFrom stats setNames
#' @importFrom SummarizedExperiment assay assayNames colData colData<-
#'
NULL

################################################################################
#' @rdname scGeneExp
#' @export
#'
scGeneExp.default <- function(scObj, gene, slot = NULL)
    stop('Unrecognized input type: scObj must be a Seurat, ',
         'SingleCellExpression, matrix or dgCMatrix object')

#' @param slot Slot where expression data is accessed. Default is 'data' for
#' Seurat objects and 'logcounts' for SingleCellExpression objects.
#'
#' @rdname scGeneExp
#' @export
#'
scGeneExp.Seurat <- function(scObj, gene, slot = 'data'){
    if (!slot %in% Layers(scObj))
        stop('Layer `', slot, '` not found for the current assay.')
    return(LayerData(scObj, slot)[gene, ])
}

#' @param slot Expression data type. Must be one of 'counts' and
#' 'logcounts'.
#'
#' @rdname scGeneExp
#' @export
#'
scGeneExp.SingleCellExperiment <- function(scObj,
                                           gene,
                                           slot = 'logcounts'){
    if (!slot %in% assayNames(scObj))
        stop('Assay `', slot, '` not found.')
    return(assay(scObj, slot)[gene, ])
}

#' @rdname scGeneExp
#' @export
#'
scGeneExp.dgCMatrix <- function(scObj, gene, slot = NULL)
    return(scObj[gene, ])

#' @rdname scGeneExp
#' @export
#'
scGeneExp.matrix <- function(scObj, gene, slot = NULL)
    return(scObj[gene, ])

################################################################################
#' @rdname scExpMat
#' @export
#'
scExpMat.default <- function(scObj,
                             slot = NULL,
                             genes = NULL,
                             densify = TRUE)
    stop('Unrecognized input type: `scObj` must be a Seurat, ',
         'SingleCellExpression, matrix or dgCMatrix object.')

#' @rdname scExpMat
#' @export
#'
scExpMat.Seurat <- function(scObj,
                            slot = 'data',
                            genes = NULL,
                            densify = TRUE){
    if (!slot %in% Layers(scObj))
        stop('Layer `', slot, '` not found for the current assay.')
    mat <- LayerData(scObj, layer=slot)
    if (!is.null(genes))
        mat <- mat[genes, ]
    if(densify)
        mat <- as.matrix(mat)
    return(mat)
}

#' @rdname scExpMat
#' @export
#'
scExpMat.SingleCellExperiment <- function(scObj,
                                          slot = 'logcounts',
                                          genes = NULL,
                                          densify = TRUE){

    if (!slot %in% assayNames(scObj))
        stop('Assay `', slot, '` not found.')
    mat <- assay(scObj, slot)
    if (!is.null(genes))
        mat <- mat[genes, ]
    if(densify)
        mat <- as.matrix(mat)
    return(mat)
}

#' @rdname scExpMat
#' @export
#'
scExpMat.dgCMatrix <- function(scObj,
                               slot = NULL,
                               genes = NULL,
                               densify = TRUE){

    if (!is.null(genes))
        mat <- scObj[genes, ]
    if(densify)
        mat <- as.matrix(mat)
    return(mat)
}

#' @rdname scExpMat
#' @export
#'
scExpMat.matrix <- function(scObj,
                            slot = NULL,
                            genes = NULL,
                            densify = TRUE){

    if (!is.null(genes))
        scObj <- scObj[genes, ]
    return(scObj)
}
