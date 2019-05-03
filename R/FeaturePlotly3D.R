#' FeaturePlotly3D
#'
#' Create a scatterplot of a given dimensional reduction set for a Seurat object,
#' coloring and sizing points by the expression level of the chosen feature.
#' Requrires a Seurat object with the reduction to be used in the corresponding
#' object@@dr slot
#'
#' @param object Seurat object
#' @param feature_use Variable to display. Currently only works with gene names
#' @param reduction_use Dimensional reduction to display. Default: umap
#' @param dim_1 Dimension to display on the x-axis. Default: 1
#' @param dim_2 Dimension to display on the y-axis. Default: 2
#' @param dim_3 Dimension to display on the z-axis. Default: 3
#' @param pt_scale Factor by which to multiply the size of the points. Default: 5
#' @param pt_shape Shape to use for the points. Default = circle
#' @param opacity Transparency level to use for the points, on a 0-1 scale. Default: 1
#' @param colors_use Color palette to use.  Accepts palettes available in the paletteer package.. Default: Reds
#' @param bins Number of bins to use in dividing expression levels.. Default: 10
#' @param plot_height Plot height in pixels. Default: 900
#' @param plot_width Plot width in pixels. Default: 900
#' @param plot_title  Display title with the name of the feature?. Default TRUE
#' @param plot_axes Display the major x, y, and z axes?. Default: FALSE
#' @param plot_grid Display the major unit tick marks?. Default: FALSE
#' @param pt_info Meta.data columns to add to the hoverinfo popup. Default: ident
#' @param legend Display legend?. Default: TRUE
#' @param legend_font_size Legend font size. Default: 12
#' @param return Return the plot object instead of displaying it. Default: FALSE
#'
#' @importFrom plotly plot_ly layout add_markers
#'
#' @return
#' @export
#'
#' @examples
FeaturePlotly3D <- function(object,
                            feature = NULL,
                            return = FALSE,
                            pt_scale = 0.5,
                            pt_shape = "circle",
                            opacity = 1,
                            reduction_use = "umap",
                            dim_1 = 1,
                            dim_2 = 2,
                            dim_3 = 3,
                            colors_use = "Red",
                            bins = 10,
                            plot_height = 900,
                            plot_width = 900,
                            plot_title = FALSE,
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

  df <- PrepInfo(object = object,
                 pt_info = pt_info,
                 df = df)

  if (is.null(feature)){
    stop("No gene or feature given")
    }

  if (is.null(feature)){ stop("No gene or feature given") }

  df <- GetFeatureValues(object = object,
                         df = df,
                         feature = feature,
                         bins = bins,
                         use.scaled = TRUE,
                         assay_use = assay_use)
  df <- GetFeatureValues(object = object,
                         df = df,
                         feature = feature,
                         use.scaled = FALSE,
                         assay_use = assay_use,
                         suffix = "size")
  df[,ncol(df)] <- df[,ncol(df)] * pt_scale

  md <- GetFeatureValues(object = object,
                         features = c(pt_info, "ident")) %>%
    mutate_at(vars(-cell),
              list(~paste0('</br> ', substitute(.), ": ", .))) %>%
    unite(info, -cell)

  df %<>% inner_join(md)

  if (isTRUE(plot_title)){
    plot_title = feature
  } else {
    plot_title = NULL
  }

  pal <- PrepQualitativePalette(bins = binds,
                                palette_use = colors_use)

  if (isTRUE(plot_title)){
    plot_title = feature
  } else {
    plot_title = NULL
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
               marker = list(symbol = pt_shape,
                             opacity = opacity,
                             sizemode = "diameter"
               ),
               width = plot_width,
               height = plot_height,
               showlegend = legend) %>%
    layout(title = plot_title,
           scene = list(
             xaxis = list(title = dim.axes[as.numeric(dim_1)], showgrid = plot_grid, visible = plot_axes),
             yaxis = list(title = dim.axes[as.numeric(dim_2)], showgrid = plot_grid, visible = plot_axes),
             zaxis = list(title = dim.axes[as.numeric(dim_3)], showgrid = plot_grid, visible = plot_axes)
           ),
           margin = c(100,NA,NA,NA)
    )

  if(!is.null(pt_info)){
    p <- p %>% add_markers(hoverinfo = "text",
                           hovertext = paste(~meta.info, ~feature),
                           showlegend = FALSE,

    )
  }

  p <- p %>% layout(legend = list(
    font = list(
      size = legend_font_size)))

  if (isTRUE(return)){
    return(df)
  } else {
    return(p)
  }
}
