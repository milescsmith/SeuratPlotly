#' FeaturePlotly3D
#'
#' Create a scatterplot of a given dimensional reduction set for a Seurat object,
#' coloring and sizing points by the expression level of the chosen feature.
#' Requrires a Seurat object with the reduction to be used in the corresponding
#' object@@dr slot
#'
#' @param object Seurat object
#' @param feature Variable to display. Currently only works with gene names
#' @param reduction Dimensional reduction to display (default: tsne)
#' @param dim_1 Dimension to display on the x-axis (default: 1)
#' @param dim_2 Dimension to display on the y-axis (default: 2)
#' @param dim_3 Dimension to display on the z-axis (default: 3)
#' @param pt_scale Factor by which to multiply the size of the points (default: 5)
#' @param pt_shape Shape to use for the points (default = circle)
#' @param opacity Transparency level to use for the points, on a 0-1 scale (default: 1)
#' @param colors Color palette to use.  Palettes from RColorBrewer and viridis. (default: Reds)
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
FeaturePlotly3D <- function(object,
                            feature = NULL,
                            return = FALSE,
                            pt_scale = 0.5,
                            pt_shape = "circle",
                            opacity = 1,
                            reduction = "tsne",
                            dim_1 = 1,
                            dim_2 = 2,
                            dim_3 = 3,
                            colors = "Reds",
                            bins = 10,
                            plot_height = 900,
                            plot_width = 900,
                            plot_title = FALSE,
                            pt_info = NULL,
                            legend = TRUE,
                            legend_font_size = 12,
                            plot_grid = FALSE,
                            plot_axes = FALSE){

  df <- as.data.frame(GetDimReduction(object = object,
                                      reduction.type = reduction,
                                      slot = "cell.embeddings"))
  dim.code <- GetDimReduction(
    object = object,
    reduction.type = reduction,
    slot = "key"
  )

  dim.axes <- colnames(
    GetDimReduction(
      object = object,
      reduction.type = reduction,
      slot = "cell.embeddings"
    )
  )

  dim.code <- c(dim.axes[[dim_1]], dim.axes[[dim_2]], dim.axes[[dim_3]])
  cell_names <- rownames(df)

  rownames(df) <- cell_names

  df$x <- df[,1]
  df$y <- df[,2]
  df$z <- df[,3]

  if(!is.null(pt_info)){
    meta.info <- list()
    # for each row
    for(i in seq(dim(df)[1])){
      # for each member of pt_info
      rowinfo = ""
      for(j in 1:length(pt_info)){
        rowinfo <- paste0(rowinfo, " </br> ", pt_info[j], ": ", object@meta.data[i, pt_info[j]])
      }
      meta.info <- c(meta.info, rowinfo)
    }
    meta.info <- unlist(meta.info)
    df$meta.info <- object@ident
  }

  if (is.null(feature)){ stop("No gene or feature given") }

  feature.data <- FetchData(object = object,
                            vars.all = feature,
                            use.scaled = TRUE)
  size.data <- FetchData(object = object,
                         vars.all = feature,
                         use.scaled = FALSE)
  feature.data[,1][feature.data[,1] == 0] <- NA
  feature.data <- as.matrix(feature.data)
  cut.feature.data <- as.numeric(as.factor(x = cut(x = as.numeric(x = feature.data), breaks = bins)))
  df[,"feature"] <- cut.feature.data
  df[,"size"] <- size.data[,1] * pt_scale

  viridis_palettes = c("viridis","inferno","magma","plasma","cividis")

  if (colors %in% rownames(brewer.pal.info)){
    pal <- colorRampPalette(brewer.pal(brewer.pal.info[colors,]$maxcolors,colors))(bins)
  } else if (colors %in% viridis_palettes){
    pal <- viridis(n = bins, option = colors)
  } else {
    pal <- colors
  }

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
    return(p)
  } else {
    p
  }
}
