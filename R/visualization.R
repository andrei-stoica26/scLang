#' @importFrom ggplot2 aes element_text ggplot geom_jitter geom_point
#' @importFrom ggplot2 geom_violin labs scale_color_gradientn
#' @importFrom ggplot2 scale_color_manual scale_fill_manual theme theme_classic
#' @importFrom henna centerTitle
#' @importFrom methods is
#' @importFrom paletteer paletteer_c
#' @importFrom wesanderson wes_palette
#'
NULL

#' Creates a dimensionality reduction plot
#'
#' This function creates a dimensionality reduction plot
#'
#' @inheritParams documentFun
#' @param noGroupsLegendLab Legend label to be used when no grouping is
#' provided (\code{groupBy} is \code{NULL})
#'
#' @return A dimensionality reduction plot.
#'
#' @examples
#' scePath <- system.file('extdata', 'sceObj.qs2', package='scLang')
#' sceObj <- qs2::qs_read(scePath)
#' dimPlot(sceObj, groupBy='Mutation_Status')
#'
#' @export
#'
dimPlot <- function(scObj,
                    groupBy = NULL,
                    title = NULL,
                    dimred = 'umap',
                    dims = c(1, 2),
                    legendTitle = 'Group',
                    noGroupsLegendLab = 'Object',
                    palette = 'grDevices::rainbow',
                    pointSize = 0.5,
                    alpha = 0.7,
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
                                 size=pointSize,
                                 alpha=alpha) +
        labs(x=paste0(dimred, '_', dims[1]),
             y=paste0(dimred, '_', dims[2]),
             color=legendTitle) +
        scale_color_manual(values=paletteer_c(palette, nColors)) +
        theme_classic() +
        theme(legend.position=legendPos,
              legend.text=element_text(size=legendTextSize),
              legend.title=element_text(size=legendTitleSize),
              axis.text=element_text(size=axisTextSize),
              axis.title=element_text(size=axisTitleSize))

    p <- centerTitle(p, title, ...)
    return(p)
}

#' Create a dimensionality reduction plot to represent a feature
#'
#' This function creates a dimensionality reduction plot to represent a
#' feature (gene expression or numeric metadata column).
#'
#' @inheritParams documentFun
#' @inheritParams pickFeature
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
                        legendTitle = NULL,
                        palette = wes_palette('Royal1')[c(3, 2)],
                        alpha = 0.6,
                        legendPos = c('right', 'top', 'left', 'bottom'),
                        legendTextSize = 10,
                        legendTitleSize = 10,
                        axisTextSize = 12,
                        axisTitleSize = 12,
                        ...){

    legendPos <- match.arg(legendPos, c('right', 'top', 'left', 'bottom'))

    if (length(dims) != 2)
        stop('`dims` must be a vector of size 2.')

    dimred <- dimredName(scObj, dimred)
    df <- as.data.frame(scDimredMat(scObj, dimred)[, dims])
    df[, 3] <- pickFeature(scObj, feature)
    colnames(df)[3] <- feature

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

#' Create a violin plot to represent a feature
#'
#' This function creates a violin plot to represent a
#' feature (gene expression or numeric metadata column).
#'
#' @inheritParams documentFun
#' @inheritParams pickFeature
#' @param xLabAngle x axis label angle.
#' @param xLabVjust x axis label vertical justification in [0, 1].
#'
#' @return A violin plot.
#'
#' @examples
#' scePath <- system.file('extdata', 'sceObj.qs2', package='scLang')
#' sceObj <- qs2::qs_read(scePath)
#' violinPlot(sceObj, 'Gene_0289')
#'
#' @export
#'
violinPlot <- function(scObj,
                       feature = rownames(scObj)[1],
                       groupBy = metadataNames(scObj)[1],
                       title = feature,
                       legendTitle = NULL,
                       xLab = 'Identity',
                       yLab = 'Expression level',
                       palette = 'grDevices::rainbow',
                       pointSize = 0.5,
                       alpha = 0.8,
                       legendPos = c('right', 'top', 'left', 'bottom'),
                       legendTextSize = 10,
                       legendTitleSize = 10,
                       axisTextSize = 12,
                       axisTitleSize = 12,
                       xLabAngle = 45,
                       xLabVjust = 0.5,
                       ...){

    legendPos <- match.arg(legendPos, c('right', 'top', 'left', 'bottom'))
    df <- data.frame(x = scCol(scObj, groupBy),
                     y = pickFeature(scObj, feature))
    nColors <- length(unique(df[, 1]))

    p <- ggplot(df, aes(x=df[, 1], y=df[, 2], fill=df[, 1])) +
        geom_violin(alpha=alpha) +
        geom_jitter(height=0, width=0.3, size=pointSize) +
        labs(x=xLab,
             y=yLab,
             fill=legendTitle) +
        scale_fill_manual(values=paletteer_c(palette, nColors)) +
        theme_classic() +
        theme(legend.position=legendPos,
              legend.text=element_text(size=legendTextSize),
              legend.title=element_text(size=legendTitleSize),
              axis.text.x=element_text(size=axisTextSize,
                                       angle=xLabAngle,
                                       vjust=xLabVjust),
              axis.text.y=element_text(size=axisTextSize),
              axis.title=element_text(size=axisTitleSize))

    p <- centerTitle(p, title, ...)
    return(p)
}
