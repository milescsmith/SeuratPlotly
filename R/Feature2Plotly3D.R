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
#' @param reduction_use Dimensional reduction to display. Default: "umap"
#' @param dim_1 Dimension to display on the x-axis. Default: 1
#' @param dim_2 Dimension to display on the y-axis. Default: 2
#' @param dim_3 Dimension to display on the z-axis. Default: 3
#' @param pt_scale Factor by which to multiply the size of the points. Default: 5
#' @param pt_shape Shape to use for the points. Default = circle
#' @param opacity Transparency level to use for the points, on a 0-1 scale. Default: 1
#' @param colors_use_1 Color palette to use for feature 1.  Palettes from RColorBrewer and viridis.. Default: Reds
#' @param colors_use_2 Color palette to use for feature 2.  Palettes from RColorBrewer and viridis.. Default: Reds
#' @param bins Number of bins to use in dividing expression levels.. Default: 10
#' @param plot_height Plot height in pixels. Default: 900
#' @param plot_width Plot width in pixels. Default: 900
#' @param plot_axes Display the major x, y, and z axes?. Default: FALSE
#' @param plot_grid Display the major unit tick marks?. Default: FALSE
#' @param pt_info Meta.data columns to add to the hoverinfo popup.. Default: ident
#' @param legend Display legend?. Default: TRUE
#' @param legend_font_size Legend font size. Default: 12
#' @param return Return the plot object instead of displaying it. Default: FALSE
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
                             reduction = "tsne",
                             dim_1 = 1,
                             dim_2 = 2,
                             dim_3 = 3,
                             colors_1 = "Reds",
                             colors_2 = "Blues",
                             bins = 10,
                             plot_height = "750",
                             plot_width = "750",
                             pt_info = NULL,
                             legend = TRUE,
                             legend_font_size = 12,
                             plot_grid = FALSE,
                             plot_axes = FALSE){

  df <- PrepDr(object,
               reduction,
               dim_1 = dim_1,
               dim_2 = dim_2,
               dim_3 = dim_3)

  #if (is.null(feature)){ stop("No gene or feature given") }

  pal1 <- PrepQuantitativePalette(bins, colors_use_1)
  pal2 <- PrepQuantitativePalette(bins, colors_use_2)

  df <- GetFeatureValues(object = object,
                         df = df,
                         feature = feature_1,
                         bins = bins,
                         use.scaled = TRUE,
                         assay_use = assay_1)
  df <- GetFeatureValues(object = object,
                         df = df,
                         feature = feature_1,
                         bins = bins,
                         use.scaled = FALSE,
                         assay_use = assay_1,
                         suffix = "size")
  df[,ncol(df)] <- df[,ncol(df)] * pt_scale

  df <- GetFeatureValues(object = object,
                         df = df,
                         feature = feature_2,
                         bins = bins,
                         use.scaled = TRUE,
                         assay_use = assay_2)
  df <- GetFeatureValues(object = object,
                         df = df,
                         feature = feature_2,
                         bins = bins,
                         use.scaled = FALSE,
                         assay_use = assay_2,
                         suffix = "size")
  df[,ncol(df)] <- df[,ncol(df)] * pt_scale

  md <- GetFeatureValues(object = object,
                         features = c(pt_info, "ident")) %>%
    mutate_at(vars(-cell),
              list(~paste0('</br> ', substitute(.), ": ", .))) %>%
    unite(info, -cell)

  df %<>% inner_join(md)

  pal_1 <- PrepQuantitativePalette(bins, colors_use_1)
  pal_2 <- PrepQuantitativePalette(bins, colors_use_2)

  p <- plot_ly(df,
               x = ~x,
               y = ~y,
               z = ~z,
               color = ~feature_1,
               mode = 'markers',
               colors = c(pal_1,pal_2),
               size = ~get(str_glue("{feature_1}_size")),
               sizes = c(0,max(df[[str_glue("{feature_1}_size")]])),
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
              color = ~feature_2,
              mode = 'markers',
              type = 'scatter3d',
              size = ~get(str_glue("{feature_2}_size")),
              sizes = c(0,max(df[[str_glue("{feature_2}_size")]])),
              marker = list(symbol = pt_shape,
                            opacity = opacity,
                            sizemode = "diameter"),
              showlegend = legend) %>%
    layout(title = str_glue("{feature_1} x {feature_2}"),
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
