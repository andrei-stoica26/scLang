#' @importFrom rlang .data
#' @importFrom S4Vectors DataFrame
#' @importFrom SeuratObject Embeddings LayerData Reductions
#' @importFrom SingleCellExperiment reducedDim reducedDims
#' @importFrom stats setNames
#' @importFrom SummarizedExperiment assay colData colData<-
#'
NULL

################################################################################
#' @rdname scGeneExp
#' @export
#'
scGeneExp.default <- function(scObj, gene, dataType = c('counts',
                                                        'data',
                                                        'logcounts'))
    stop('Unrecognized input type: scObj must be a Seurat, ',
         'SingleCellExpression, matrix or dgCMatrix object')

#' @param dataType Expression data type. Must be one of 'counts' and 'data'.
#'
#' @rdname scGeneExp
#' @export
#'
scGeneExp.Seurat <- function(scObj,
                             gene,
                             dataType = c('counts',
                                          'data',
                                          'logcounts')){
    dataType <- match.arg(dataType, c('counts', 'data', 'logcounts'))
    if (dataType == 'logcounts')
        dataType <- 'data'
    return(LayerData(scObj, dataType)[gene, ])
}

#' @param dataType Expression data type. Must be one of 'counts' and
#' 'logcounts'.
#'
#' @rdname scGeneExp
#' @export
#'
scGeneExp.SingleCellExperiment <- function(scObj,
                                           gene,
                                           dataType = c('counts',
                                                        'data',
                                                        'logcounts')){
    dataType <- match.arg(dataType, c('counts', 'data', 'logcounts'))
    if (dataType == 'data')
        dataType <- 'logcounts'
    return(assay(scObj, dataType)[gene, ])
}

#' @rdname scGeneExp
#' @export
#'
scGeneExp.dgCMatrix <- function(scObj,
                                gene,
                                dataType = c('counts',
                                             'data',
                                             'logcounts'))
    return(scObj[gene, ])

#' @rdname scGeneExp
#' @export
#'
scGeneExp.matrix <- function(scObj,
                             gene,
                             dataType = c('counts',
                                          'data',
                                          'logcounts'))
    return(scObj[gene, ])

################################################################################
#' @rdname scExpMat
#' @export
#'
scExpMat.default <- function(scObj,
                             dataType = c('data',
                                          'counts',
                                          'logcounts'),
                             genes = NULL,
                             densify = TRUE)
    stop('Unrecognized input type: `scObj` must be a Seurat, ',
         'SingleCellExpression, matrix or dgCMatrix object.')

#' @rdname scExpMat
#' @export
#'
scExpMat.Seurat <- function(scObj,
                            dataType = c('data',
                                         'counts',
                                         'logcounts'),
                            genes = NULL,
                            densify = TRUE){
    dataType <- match.arg(dataType, c('data', 'counts', 'logcounts'))
    if (dataType == 'logcounts')
        dataType <- 'data'
    mat <- LayerData(scObj, layer=dataType)
    if (!is.null(genes))
        mat <- mat[genes, ]
    if(densify)
        mat <- suppressWarnings(as.matrix(mat))
    return(mat)
}

#' @rdname scExpMat
#' @export
#'
scExpMat.SingleCellExperiment <- function(scObj,
                                          dataType = c('data',
                                                       'counts',
                                                       'logcounts'),
                                          genes = NULL,
                                          densify = TRUE){

    dataType <- match.arg(dataType, c('data', 'counts', 'logcounts'))
    if (dataType == 'data')
        dataType <- 'logcounts'
    mat <- assay(scObj, dataType)
    if (!is.null(genes))
        mat <- mat[genes, ]
    if(densify)
        mat <- suppressWarnings(as.matrix(mat))
    return(mat)
}

#' @rdname scExpMat
#' @export
#'
scExpMat.dgCMatrix <- function(scObj,
                               dataType = c('data',
                                            'counts',
                                            'logcounts'),
                               genes = NULL,
                               densify = TRUE){

    if (!is.null(genes))
        mat <- scObj[genes, ]
    if(densify)
        mat <- suppressWarnings(as.matrix(mat))
    return(mat)
}

#' @rdname scExpMat
#' @export
#'
scExpMat.matrix <- function(scObj,
                            dataType = c('data',
                                         'counts',
                                         'logcounts'),
                            genes = NULL,
                            densify = TRUE){

    if (!is.null(genes))
        scObj <- scObj[genes, ]
    return(scObj)
}
