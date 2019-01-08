#' Feature2Plotly3D
#'
#' Create a scatterplot of a given dimensional reduction set for a Seurat object,
#' coloring and sizing points by the expression level of the chosen feature.
#' Requrires a Seurat object with the reduction to be used in the corresponding
#' seuratObj@@dr slot
#'
#' @param seuratObj Seurat object
#' @param feature.1.use First variable to display. Currently only works with gene names
#' @param feature.2.use Second variable to display. Currently only works with gene names
#' @param reduction.use Dimensional reduction to display (default: tsne)
#' @param dim.1 Dimension to display on the x-axis (default: 1)
#' @param dim.2 Dimension to display on the y-axis (default: 2)
#' @param dim.3 Dimension to display on the z-axis (default: 3)
#' @param pt.scale Factor by which to multiply the size of the points (default: 5)
#' @param pt.shape Shape to use for the points (default = circle)
#' @param opacity Transparency level to use for the points, on a 0-1 scale (default: 1)
#' @param colors.use.1 Color palette to use for feature 1.  Palettes from RColorBrewer and viridis. (default: Reds)
#' @param colors.use.2 Color palette to use for feature 2.  Palettes from RColorBrewer and viridis. (default: Reds)
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
Feature2Plotly3D <- function(seuratObj,
                             feature.1.use = NULL,
                             feature.2.use = NULL,
                             do.return = FALSE,
                             pt.scale = 0.5,
                             pt.shape = "circle",
                             opacity = 1,
                             reduction.use = "tsne",
                             dim.1 = 1,
                             dim.2 = 2,
                             dim.3 = 3,
                             colors.use.1 = "Reds",
                             colors.use.2 = "Blues",
                             bins = 10,
                             plot.height = "750",
                             plot.width = "750",
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
  df <- df[,dim.code]
  cell_names <- rownames(df)

  df$x <- df[,1]
  df$y <- df[,2]
  df$z <- df[,3]
  rownames(df) <- cell_names

  #if (is.null(feature)){ stop("No gene or feature given") }

  feature.1.data <- FetchData(object = seuratObj,
                              vars.all = feature.1.use,
                              use.scaled = TRUE)
  size.1.data <- FetchData(object = seuratObj,
                           vars.all = feature.1.use,
                           use.scaled = FALSE)
  feature.1.data <- as.matrix(feature.1.data)
  cut.feature.1.data <- as.numeric(as.factor(x = cut(x = as.numeric(x = feature.1.data), breaks = bins)))
  cut.feature.1.data <- as.factor(as.numeric(as.factor(x = cut(x = as.numeric(x = feature.1.data), breaks = bins))))
  df[,"feature.1"] <- cut.feature.1.data
  df[,"size.1"] <- size.1.data[,1] * pt.scale

  feature.2.data <- FetchData(object = seuratObj,
                              vars.all = feature.2.use,
                              use.scaled = TRUE)
  size.2.data <- FetchData(object = seuratObj,
                           vars.all = feature.2.use,
                           use.scaled = FALSE)
  feature.2.data[,1][feature.2.data[,1] == 0] <- NA
  feature.2.data <- as.matrix(feature.2.data)
  cut.feature.2.data <- as.numeric(as.factor(x = cut(x = as.numeric(x = feature.2.data), breaks = bins)))
  df[,"feature.2"] <- cut.feature.2.data
  df[,"size.2"] <- size.2.data[,1] * pt.scale

  if(!is.null(pt.info)){
    meta.info <- list()
    # for each row
    for(i in seq(dim(df)[1])){
      # for each member of pt.info
      rowinfo = ""
      for(j in 1:length(pt.info)){
        rowinfo <- paste0(rowinfo, " </br> ", pt.info[j], ": ", seuratObj@meta.data[i, pt.info[j]])
      }
      rowinfo <- paste0(rowinfo, ' </br> Expr ', feature.1.use, ': ', df[i,"feature.1"])
      rowinfo <- paste0(rowinfo, ' </br> Expr ', feature.2.use, ': ', df[i,"feature.2"])
      meta.info <- c(meta.info, rowinfo)
    }
    meta.info <- unlist(meta.info)
    df$meta.info <- meta.info
  }

  viridis_palettes = c("viridis","inferno","magma","plasma","cividis")

  if (colors.use.1 %in% rownames(brewer.pal.info)){
    pal.1 <- colorRampPalette(brewer.pal(brewer.pal.info[colors.use.1,]$maxcolors,colors.use.1))(bins)
  } else if (colors.use.1 %in% viridis_palettes){
    pal.1 <- viridis(n = bins, option = colors.use.1)
  } else {
    pal.1 <- colors.use.1
  }

  if (colors.use.2 %in% rownames(brewer.pal.info)){
    pal.2 <- colorRampPalette(brewer.pal(brewer.pal.info[colors.use.2,]$maxcolors,colors.use.2))(bins)
  } else if (colors.use.2 %in% viridis_palettes){
    pal.2 <- viridis(n = bins, option = colors.use.2)
  } else {
    pal.2 <- colors.use.2
  }

  p <- plot_ly(df,
               x = ~x,
               y = ~y,
               z = ~z,
               color = ~feature.1,
               mode = 'markers',
               colors = c(pal.1,pal.2),
               size = ~size.1,
               sizes = c(0,max(df$size.1)),
               marker = list(symbol = pt.shape,
                             opacity = opacity,
                             sizemode = "diameter"),
               width = plot.width,
               height = plot.height,
               type = 'scatter3d',
               showlegend = legend) %>%
    add_trace(df,
              x = ~x,
              y = ~y,
              z = ~z,
              color = ~feature.2,
              mode = 'markers',
              type = 'scatter3d',
              size = ~size.2,
              sizes = c(0,max(df$size.2)),
              marker = list(symbol = pt.shape,
                            opacity = opacity,
                            sizemode = "diameter"),
              showlegend = legend) %>%
    layout(title = paste(feature.1.use, feature.2.use, sep = " x "),
           scene = list(
             aspectratio = list(x = 1,y = 1,z = 1),
             camera = list(
               center = list(x = 0,y = 0,z = 0),
               eye = list(x = 2,y = -1,z = 0.5),
               up = list(x = 0,y = 0,z = 1)),
             dragmode = "turnable",
             xaxis = list(title = dim.axes[as.numeric(dim.1)], type = "double", showgrid = plot.grid, visible = plot.axes),
             yaxis = list(title = dim.axes[as.numeric(dim.2)], type = "double", showgrid = plot.grid, visible = plot.axes),
             zaxis = list(title = dim.axes[as.numeric(dim.3)], type = "double", showgrid = plot.grid, visible = plot.axes),
             margin = c(100,NA,NA,NA)))

  if(!is.null(pt.info)){
    p <- p %>% add_markers(hoverinfo = "text",
                           hovertext = ~meta.info,
                           showlegend = FALSE
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
