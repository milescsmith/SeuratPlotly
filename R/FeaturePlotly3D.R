#' FeaturePlotly3D
#'
#' Create a scatterplot of a given dimensional reduction set for a Seurat object,
#' coloring and sizing points by the expression level of the chosen feature.
#' Requrires a Seurat object with the reduction to be used in the corresponding
#' seuratObj@@dr slot
#'
#' @param seuratObj Seurat object
#' @param feature.use Variable to display. Currently only works with gene names
#' @param reduction.use Dimensional reduction to display (default: tsne)
#' @param dim.1 Dimension to display on the x-axis (default: 1)
#' @param dim.2 Dimension to display on the y-axis (default: 2)
#' @param dim.3 Dimension to display on the z-axis (default: 3)
#' @param pt.scale Factor by which to multiply the size of the points (default: 5)
#' @param pt.shape Shape to use for the points (default = circle)
#' @param opacity Transparency level to use for the points, on a 0-1 scale (default: 1)
#' @param colors.use Color palette to use.  Palettes from RColorBrewer and viridis. (default: Reds)
#' @param bins Number of bins to use in dividing expression levels. (default: 10)
#' @param plot.height Plot height in pixels (default: 900)
#' @param plot.width Plot width in pixels (default: 900)
#' @param plot.title  Display title with the name of the feature? (default TRUE)
#' @param plot.axes Display the major x, y, and z axes? (default: FALSE)
#' @param plot.grid Display the major unit tick marks? (default: FALSE)
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
#' @return If do.return is TRUE, a plotly object.
#' @export
#'
#' @examples
FeaturePlotly3D <- function(seuratObj,
                            feature.use = NULL,
                            do.return = FALSE,
                            pt.scale = 0.5,
                            pt.shape = "circle",
                            opacity = 1,
                            reduction.use = "tsne",
                            dim.1 = 1,
                            dim.2 = 2,
                            dim.3 = 3,
                            colors.use = "Reds",
                            bins = 10,
                            plot.height = 900,
                            plot.width = 900,
                            plot.title = FALSE,
                            pt.info = NULL,
                            legend = TRUE,
                            legend.font.size = 12,
                            plot.grid = FALSE,
                            plot.axes = FALSE){

  df <- as.data.frame(GetDimReduction(object = seuratObj,
                                      reduction.type = reduction.use,
                                      slot = "cell.embeddings"))
  dim.code <- GetDimReduction(
    object = seuratObj,
    reduction.type = reduction.use,
    slot = "key"
  )

  dim.axes <- colnames(
    GetDimReduction(
      object = seuratObj,
      reduction.type = reduction.use,
      slot = "cell.embeddings"
    )
  )

  dim.code <- c(dim.axes[[dim.1]], dim.axes[[dim.2]], dim.axes[[dim.3]])
  cell_names <- rownames(df)

  rownames(df) <- cell_names

  df$x <- df[,1]
  df$y <- df[,2]
  df$z <- df[,3]

  if(!is.null(pt.info)){
    meta.info <- list()
    # for each row
    for(i in seq(dim(df)[1])){
      # for each member of pt.info
      rowinfo = ""
      for(j in 1:length(pt.info)){
        rowinfo <- paste0(rowinfo, " </br> ", pt.info[j], ": ", seuratObj@meta.data[i, pt.info[j]])
      }
      meta.info <- c(meta.info, rowinfo)
    }
    meta.info <- unlist(meta.info)
    df$meta.info <- seuratObj@ident
  }

  if (is.null(feature.use)){ stop("No gene or feature given") }

  feature.data <- FetchData(object = seuratObj,
                            vars.all = feature.use,
                            use.scaled = TRUE)
  size.data <- FetchData(object = seuratObj,
                         vars.all = feature.use,
                         use.scaled = FALSE)
  feature.data[,1][feature.data[,1] == 0] <- NA
  feature.data <- as.matrix(feature.data)
  cut.feature.data <- as.numeric(as.factor(x = cut(x = as.numeric(x = feature.data), breaks = bins)))
  df[,"feature"] <- cut.feature.data
  df[,"size"] <- size.data[,1] * pt.scale

  viridis_palettes = c("viridis","inferno","magma","plasma","cividis")

  if (colors.use %in% rownames(brewer.pal.info)){
    pal <- colorRampPalette(brewer.pal(brewer.pal.info[colors.use,]$maxcolors,colors.use))(bins)
  } else if (colors.use %in% viridis_palettes){
    pal <- viridis(n = bins, option = colors.use)
  } else {
    pal <- colors.use
  }

  if (isTRUE(plot.title)){
    plot.title = feature.use
  } else {
    plot.title = NULL
  }

  p <- plot_ly(df,
               x = ~x,
               y = ~y,
               z = ~z,
               color = ~feature,
               mode = 'markers',
               type = "scatter3d",
               colors = pal,
               size = ~size,
               sizes = c(0,max(df$size)),
               marker = list(symbol = pt.shape,
                             opacity = opacity,
                             sizemode = "diameter"
               ),
               width = plot.width,
               height = plot.height,
               showlegend = legend) %>%
    layout(title = plot.title,
           scene = list(
             xaxis = list(title = dim.axes[as.numeric(dim.1)], showgrid = plot.grid, visible = plot.axes),
             yaxis = list(title = dim.axes[as.numeric(dim.2)], showgrid = plot.grid, visible = plot.axes),
             zaxis = list(title = dim.axes[as.numeric(dim.3)], showgrid = plot.grid, visible = plot.axes)
           ),
           margin = c(100,NA,NA,NA)
    )

  if(!is.null(pt.info)){
    p <- p %>% add_markers(hoverinfo = "text",
                           hovertext = paste(~meta.info, ~feature),
                           showlegend = FALSE,

    )
  }

  p <- p %>% layout(legend = list(
    font = list(
      size = legend.font.size)))

  if (isTRUE(do.return)){
    return(p)
  } else {
    p
  }
}
