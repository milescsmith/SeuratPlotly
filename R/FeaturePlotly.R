#' FeaturePlotly
#'
#' Create a scatterplot of a given dimensional reduction set for a Seurat object,
#' coloring and sizing points by the expression level of the chosen feature.
#' Requrires a Seurat object with the reduction to be used in the corresponding
#' object@@dr slot
#'
#' @param object Seurat object
#' @param feature Variable to display. Currently only works with gene names
#' @param reduction.use Dimensional reduction to display (default: tsne)
#' @param dim.1 Dimension to display on the x-axis (default: 1)
#' @param dim.2 Dimension to display on the y-axis (default: 2)
#' @param pt.scale Factor by which to multiply the size of the points (default: 5)
#' @param pt.shape Shape to use for the points (default = circle)
#' @param opacity Transparency level to use for the points, on a 0-1 scale (default: 1)
#' @param colors.use Color palette to use.  Palettes from RColorBrewer and viridis or a list of colors. (default: Reds)
#' @param bins Number of bins to use in dividing expression levels. (default: 10)
#' @param plot.height Plot height in pixels (default: 900)
#' @param plot.width Plot width in pixels (default: 900)
#' @param plot.title  Display title with the name of the feature? (default TRUE)
#' @param pt.info Meta.data columns to add to the hoverinfo popup. (default: ident)
#' @param legend Display legend? (default: TRUE)
#' @param legend.font.size Legend font size (default: 12)
#' @param do.return Return the plot object instead of displaying it (default: FALSE)
#'
#' @import dplyr
#' @importFrom magrittr "%>%"
#' @importFrom Seurat GetDimReduction
#' @importFrom Seurat FetchData
#' @importFrom RColorBrewer brewer.pal
#' @importFrom RColorBrewer brewer.pal.info
#' @importFrom viridis viridis
#' @importFrom plotly plot_ly
#' @importFrom plotly layout
#' @importFrom grDevices colorRampPalette
#'
#' @return plotly object
#' @export
#'
#' @examples
FeaturePlotly <- function(object,
                          feature = NULL,
                          assay.use = "RNA",
                          slot.use = "data",
                          do.return = TRUE,
                          pt.scale = 5,
                          pt.shape = "circle",
                          opacity = 1,
                          reduction.use = "tsne",
                          dim.1 = 1,
                          dim.2 = 2,
                          colors.use = c("blue","red"),
                          reverse.color.scale = FALSE,
                          bins = 10,
                          plot.height = 900,
                          plot.width = 900,
                          plot.title = FALSE,
                          pt.info = NULL,
                          legend = TRUE,
                          legend.font.size = 12){

  df <- PrepDf(object,
               reduction.use,
               dim.1 = dim.1,
               dim.2 = dim.2)

  df <- PrepInfo(object = object,
                 pt.info = c(feature),
                 df = df)

  if (is.null(feature)){ stop("No gene or feature given") }

  df <- GetFeatureValues(object = object,
                         df = df,
                         feature = feature,
                         assay.use,
                         slot.use)

  cut.feature.data <- as.numeric(as.factor(x = cut(x = as.numeric(df[,feature]), breaks = bins)))

  df[,"size"] <- df[,feature] * pt.scale

  if (isTRUE(plot.title)){
    plot.title = feature
  } else {
    plot.title = NULL
  }

  p <- plot_ly(df,
               x = ~x,
               y = ~y,
               mode = 'markers',
               type = "scattergl",
               size = ~size,
               sizes = c(0,max(df$size)),
               marker = list(symbol = pt.shape,
                             opacity = opacity,
                             color = ~get(feature),
                             line = list(width = 0),
                             colorscale=colors.use,
                             reversescale = reverse.color.scale,
                             cmin = 0,
                             cmax = 1,
                             sizemode = "diameter",
                             colorbar = list(title = feature)),
               width = plot.width,
               height = plot.height,
               showlegend = legend,
               hoverinfo = "text",
               text = ~meta.info) %>%
    layout(
      title = plot.title,
      xaxis = list(title = glue("{reduction.use}_{dim.1}")),
      yaxis = list(title = glue("{reduction.use}_{dim.2}"))
    )

  p <- p %>% layout(legend = list(
    font = list(
      size = legend.font.size)
    )
  )

  if (isTRUE(do.return)){
    return(p)
  } else {
    p
  }
}
