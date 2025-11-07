#' @importFrom ggplot2 ggplot geom_point
#' @importFrom henna centerTitle
#'
NULL

#' Creates a dimensionality reduction plot
#'
#' This function creates a dimensionality reduction plot
#'
#' @param title Plot title.
#' @param paletteFun Palette function. Must accept the number of colors as the
#' first argument and require no other arguments.
#'
#' @return A dimensionality reduction plot
#'
#' @export
#'
dimPlot <- function(scObj,
                    title = NULL,
                    groupBy = NULL,
                    dimred = 'umap',
                    dims = c(1, 2),
                    legendTitle = 'Group',
                    noGroupsLegendLab = 'Object',
                    paletteFun = rainbow,
                    alpha=0.7){

    dimred <- dimredName(scObj, dimred)
    df <- as.data.frame(scDimredMat(scObj, dimred)[, dims])
    if(is.null(groupBy)){
        df[[noGroupsLegendLab]] <- noGroupsLegendLab
        legendTitle <- NULL
    } else
        df[[groupBy]] <- scCol(scObj, groupBy)

    nColors <- length(unique(df[, 3]))
    p <- ggplot(df) + geom_point(aes(x=df[, 1],
                                     y=df[, 2],
                                     color=df[, 3]),
                                 alpha=alpha) +
        labs(x=paste0(dimred, '_', dims[1]),
             y=paste0(dimred, '_', dims[2]),
             color=legendTitle) +
        scale_color_manual(values=paletteFun(nColors)) +
        theme_classic()

    p <- centerTitle(p, title)
    return(p)
}
