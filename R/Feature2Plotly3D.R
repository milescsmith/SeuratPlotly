#' Feature2Plotly3D
#'
#' Create a scatterplot of a given dimensional reduction set for a Seurat object,
#' coloring and sizing points by the expression level of the chosen feature.
#' Requrires a Seurat object with the reduction to be used in the corresponding
#' object@@dr slot
#'
#' @param object Seurat object
#' @param feature_1 First variable to display. Currently only works with gene names
#' @param feature_2 Second variable to display. Currently only works with gene names
#' @param reduction_use Dimensional reduction to display (default: tsne)
#' @param dim_1 Dimension to display on the x-axis (default: 1)
#' @param dim_2 Dimension to display on the y-axis (default: 2)
#' @param dim_3 Dimension to display on the z-axis (default: 3)
#' @param pt_scale Factor by which to multiply the size of the points (default: 5)
#' @param pt_shape Shape to use for the points (default = circle)
#' @param opacity Transparency level to use for the points, on a 0-1 scale (default: 1)
#' @param colors_use_1 Color palette to use for feature 1.  Palettes from RColorBrewer and viridis. (default: Reds)
#' @param colors_use_2 Color palette to use for feature 2.  Palettes from RColorBrewer and viridis. (default: Reds)
#' @param bins Number of bins to use in dividing expression levels. (default: 10)
#' @param plot_height Plot height in pixels (default: 900)
#' @param plot_width Plot width in pixels (default: 900)
#' @param plot_title  Display title with the name of the feature? (default TRUE)
#' @param plot_axes Display the major x, y, and z axes? (default: FALSE)
#' @param plot_grid Display the major unit tick marks? (default: FALSE)
#' @param pt_info Meta.data columns to add to the hoverinfo popup. (default: ident)
#' @param legend Display legend? (default: TRUE)
#' @param legend_font_size Legend font size (default: 12)
#' @param return Return the plot object instead of displaying it (default: FALSE)
#'
#' @import dplyr
#' @import Seurat
#' @importFrom RColorBrewer brewer.pal brewer.pal.info
#' @importFrom viridis viridis
#' @importFrom plotly plot_ly layout
#' @importFrom grDevices colorRampPalette
#'
#' @return If return is TRUE, a plotly object.
#' @export
#'
#' @examples
Feature2Plotly3D <- function(object,
                             feature_1 = NULL,
                             feature_2 = NULL,
                             return = FALSE,
                             pt_scale = 0.5,
                             pt_shape = "circle",
                             opacity = 1,
                             reduction_use = "tsne",
                             dim_1 = 1,
                             dim_2 = 2,
                             dim_3 = 3,
                             colors_use_1 = "Reds",
                             colors_use_2 = "Blues",
                             bins = 10,
                             plot_height = "750",
                             plot_width = "750",
                             pt_info = NULL,
                             legend = TRUE,
                             legend_font_size = 12,
                             plot_grid = FALSE,
                             plot_axes = FALSE){

  df <- as.data.frame(GetDimReduction(object = object,
                                      reduction.type = reduction_use,
                                      slot = "cell.embeddings"))
  dim.code <- GetDimReduction(
    object = object,
    reduction.type = reduction_use,
    slot = "key"
  )

  dim.axes <- colnames(
    GetDimReduction(
      object = object,
      reduction.type = reduction_use,
      slot = "cell.embeddings"
    )
  )

  dim.code <- c(dim.axes[[dim_1]], dim.axes[[dim_2]], dim.axes[[dim_3]])
  df <- df[,dim.code]
  cell_names <- rownames(df)

  df$x <- df[,1]
  df$y <- df[,2]
  df$z <- df[,3]
  rownames(df) <- cell_names

  #if (is.null(feature)){ stop("No gene or feature given") }

  feature_1_data <- FetchData(object = object,
                              vars.all = feature_1,
                              use.scaled = TRUE)
  size.1.data <- FetchData(object = object,
                           vars.all = feature_1,
                           use.scaled = FALSE)
  feature_1_data <- as.matrix(feature_1_data)
  cut.feature_1_data <- as.numeric(as.factor(x = cut(x = as.numeric(x = feature_1_data), breaks = bins)))
  cut.feature_1_data <- as.factor(as.numeric(as.factor(x = cut(x = as.numeric(x = feature_1_data), breaks = bins))))
  df[,"feature.1"] <- cut.feature_1_data
  df[,"size.1"] <- size.1.data[,1] * pt_scale

  feature_2_data <- FetchData(object = object,
                              vars.all = feature_2,
                              use.scaled = TRUE)
  size.2.data <- FetchData(object = object,
                           vars.all = feature_2,
                           use.scaled = FALSE)
  feature_2_data[,1][feature_2_data[,1] == 0] <- NA
  feature_2_data <- as.matrix(feature_2_data)
  cut.feature_2_data <- as.numeric(as.factor(x = cut(x = as.numeric(x = feature_2_data), breaks = bins)))
  df[,"feature.2"] <- cut.feature_2_data
  df[,"size.2"] <- size.2.data[,1] * pt_scale

  if(!is.null(pt_info)){
    meta.info <- list()
    # for each row
    for(i in seq(dim(df)[1])){
      # for each member of pt_info
      rowinfo = ""
      for(j in 1:length(pt_info)){
        rowinfo <- paste0(rowinfo, " </br> ", pt_info[j], ": ", object@meta.data[i, pt_info[j]])
      }
      rowinfo <- paste0(rowinfo, ' </br> Expr ', feature_1, ': ', df[i,"feature.1"])
      rowinfo <- paste0(rowinfo, ' </br> Expr ', feature_2, ': ', df[i,"feature.2"])
      meta.info <- c(meta.info, rowinfo)
    }
    meta.info <- unlist(meta.info)
    df$meta.info <- meta.info
  }

  viridis_palettes = c("viridis","inferno","magma","plasma","cividis")

  if (colors_use_1 %in% rownames(brewer.pal.info)){
    pal.1 <- colorRampPalette(brewer.pal(brewer.pal.info[colors_use_1,]$maxcolors,colors_use_1))(bins)
  } else if (colors_use_1 %in% viridis_palettes){
    pal.1 <- viridis(n = bins, option = colors_use_1)
  } else {
    pal.1 <- colors_use_1
  }

  if (colors_use_2 %in% rownames(brewer.pal.info)){
    pal.2 <- colorRampPalette(brewer.pal(brewer.pal.info[colors_use_2,]$maxcolors,colors_use_2))(bins)
  } else if (colors_use_2 %in% viridis_palettes){
    pal.2 <- viridis(n = bins, option = colors_use_2)
  } else {
    pal.2 <- colors_use_2
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
               marker = list(symbol = pt_shape,
                             opacity = opacity,
                             sizemode = "diameter"),
               width = plot_width,
               height = plot_height,
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
              marker = list(symbol = pt_shape,
                            opacity = opacity,
                            sizemode = "diameter"),
              showlegend = legend) %>%
    layout(title = paste(feature_1, feature_2, sep = " x "),
           scene = list(
             aspectratio = list(x = 1,y = 1,z = 1),
             camera = list(
               center = list(x = 0,y = 0,z = 0),
               eye = list(x = 2,y = -1,z = 0.5),
               up = list(x = 0,y = 0,z = 1)),
             dragmode = "turnable",
             xaxis = list(title = dim.axes[as.numeric(dim_1)], type = "double", showgrid = plot_grid, visible = plot_axes),
             yaxis = list(title = dim.axes[as.numeric(dim_2)], type = "double", showgrid = plot_grid, visible = plot_axes),
             zaxis = list(title = dim.axes[as.numeric(dim_3)], type = "double", showgrid = plot_grid, visible = plot_axes),
             margin = c(100,NA,NA,NA)))

  if(!is.null(pt_info)){
    p <- p %>% add_markers(hoverinfo = "text",
                           hovertext = ~meta.info,
                           showlegend = FALSE
    )
  }
  p <- p %>% layout(legend = list(
    font = list(
      size = legend_font_size)))

  if (isTRUE(return)){
    return(p)
  } else {
    p
  }
}
