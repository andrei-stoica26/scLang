###############################################################################
#' @rdname metadataDF
#' @export
#'
metadataDF.default <- function(scObj)
    stop('Unrecognized input type: `scObj` must be a Seurat or ',
         'SingleCellExpression object.')

#' @rdname metadataDF
#' @export
#'
`metadataDF<-.default` <- function(scObj, value)
    stop('Unrecognized input type: `scObj` must be a Seurat or ',
         'SingleCellExpression object.')

#' @rdname metadataDF
#' @export
#'
metadataDF.Seurat <- function(scObj)
    return(scObj[[]])

#' @rdname metadataDF
#' @export
#'
`metadataDF<-.Seurat` <- function(scObj, value){
    if (!is.data.frame(value))
        stop('`value` must be a data.frame.')
    scObj[[]] <- value
    return(scObj)
}

#' @rdname metadataDF
#' @export
#'
metadataDF.SingleCellExperiment <- function(scObj)
    return(as.data.frame(colData(scObj)))

#' @rdname metadataDF
#' @export
#'
`metadataDF<-.SingleCellExperiment` <- function(scObj, value) {
    if (!is.data.frame(value))
        stop('`value` must be a data.frame.')
    colData(scObj) <- DataFrame(value)
    return(scObj)
}

###############################################################################
#' Return metadata names
#'
#' This function extracts metadata names from a Seurat or
#' SingleCellExperiment object. It can also be used to modify metadata names.
#'
#' @inheritParams metadataDF
#'
#' @return The names of the metadata columns.
#'
#' @examples
#' scePath <- system.file('extdata', 'sceObj.qs2', package='scLang')
#' sceObj <- qs2::qs_read(scePath)
#' colNames <- metadataNames(sceObj)
#'
#' @export
#'
metadataNames <- function(scObj)
    return(colnames(metadataDF(scObj)))

#' @rdname metadataNames
#' @export
#'
`metadataNames<-` <- function(scObj, value){
    if (!is.character(value))
        stop('`value` must be a character.')
    colnames(metadataDF(scObj)) <- value
    return(scObj)
}

###############################################################################
#' @rdname scCol
#' @export
#'
scCol.default <- function(scObj, col)
    stop('Unrecognized input type: `scObj` must be a Seurat or ',
         'SingleCellExpression object')

#' @rdname scCol
#' @export
#'
`scCol<-.default` <- function(scObj, col, value)
    stop('Unrecognized input type: `scObj` must be a Seurat or ',
         'SingleCellExpression object.')

#' @rdname scCol
#' @export
#'
scCol.Seurat <- function(scObj, col)
    return(scObj[[]][[col]])

#' @rdname scCol
#' @export
#'
`scCol<-.Seurat` <- function(scObj, col, value){
    scObj[[]][[col]] <- value
    return(scObj)
}

#' @rdname scCol
#' @export
#'
scCol.SingleCellExperiment <- function(scObj, col)
    return(scObj[[col]])

#' @rdname scCol
#' @export
#'
`scCol<-.SingleCellExperiment` <- function(scObj, col, value){
    scObj[[col]] <- value
    return(scObj)
}
