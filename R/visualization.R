#' @importFrom ggplot2 aes element_text ggplot geom_point labs scale_color_gradientn scale_color_manual theme theme_classic
#' @importFrom grDevices rainbow
#' @importFrom henna centerTitle
#' @importFrom methods is
#' @importFrom wesanderson wes_palette
#'
NULL

#' Creates a dimensionality reduction plot
#'
#' This function creates a dimensionality reduction plot
#'
#' @inheritParams documentFun
#' @param paletteFun Palette function. Must accept the number of colors as the
#' first argument and require no other arguments.
#' @param noGroupsLegendLab Legend label to be used when no grouping is
#' provided (\code{groupBy} is \code{NULL})
#'
#' @return A dimensionality reduction plot.
#'
#' @examples
#' scePath <- system.file('extdata', 'sceObj.qs2', package='scLang')
#' sceObj <- qs2::qs_read(scePath)
#' dimPlot(sceObj, groupBy='Donor')
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
                    alpha=0.6,
                    legendPos = c('right', 'top', 'left', 'bottom'),
                    legendTextSize = 10,
                    legendTitleSize = 10,
                    axisTextSize = 12,
                    axisTitleSize = 12,
                    ...
                    ){
    legendPos <- match.arg(legendPos, c('right', 'top', 'left', 'bottom'))
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
        theme_classic() +
        theme(legend.position=legendPos,
              legend.text=element_text(size=legendTextSize),
              legend.title=element_text(size=legendTitleSize),
              axis.text=element_text(size=axisTextSize),
              axis.title=element_text(size=axisTitleSize))

    p <- centerTitle(p, title, ...)
    return(p)
}

#' Creates a dimensionality reduction plot
#'
#' This function creates a dimensionality reduction plot
#'
#' @inheritParams documentFun
#' @inheritParams pickFeature
#' @param palette Color palette.
#'
#' @return A feature plot.
#'
#' @examples
#' scePath <- system.file('extdata', 'sceObj.qs2', package='scLang')
#' sceObj <- qs2::qs_read(scePath)
#' featurePlot(sceObj, 'Gene_0289')
#'
#' @export
#'
featurePlot <- function(scObj,
                        feature = rownames(scObj)[1],
                        title = feature,
                        dimred = 'umap',
                        dims = c(1, 2),
                        legendTitle = feature,
                        palette = wes_palette('Royal1')[c(3, 2)],
                        alpha = 0.6,
                        legendPos = c('right', 'top', 'left', 'bottom'),
                        legendTextSize = 10,
                        legendTitleSize = 10,
                        axisTextSize = 12,
                        axisTitleSize = 12,
                        ...){

    legendPos <- match.arg(legendPos, c('right', 'top', 'left', 'bottom'))

    df <- pickFeature(scObj, feature, dimred, dims)

    p <- ggplot(df) + geom_point(aes(x=df[, 1],
                                     y=df[, 2],
                                     color=df[, 3]),
                                 alpha=alpha) +
        labs(x=paste0(dimred, '_', dims[1]),
             y=paste0(dimred, '_', dims[2]),
             color=legendTitle) +
        scale_color_gradientn(colors=palette) +
        theme_classic() +
        theme(legend.position=legendPos,
              legend.text=element_text(size=legendTextSize),
              legend.title=element_text(size=legendTitleSize),
              axis.text=element_text(size=axisTextSize),
              axis.title=element_text(size=axisTitleSize))

    p <- centerTitle(p, title, ...)
    return(p)
}

