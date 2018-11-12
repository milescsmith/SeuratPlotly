#' DimPlotly3D
#'
#' @param seuratObj Seurat object
#' @param group.by Variable by which to group cells. Currently only works with the current ident and column names from meta.data (default: ident)
#' @param reduction.use Dimensional reduction to display (default: tsne)
#' @param dim.1 Dimension to display on the x-axis (default: 1)
#' @param dim.2 Dimension to display on the y-axis (default: 2)
#' @param dim.3 Dimension to display on the z-axis (default: 3)
#' @param do.label Add a label showing thr group name to the graph (default: FALSE)
#' @param label.size Label font size (default: 12)
#' @param show.arrow Offset the position of the labels and instead point to each group with an arrow (default: FALSE)
#' @param label.color Color for label border and arrow.  Need hex value. (default = '000000')
#' @param do.return Return the plot object instead of displaying it (default: FALSE)
#' @param pt.size Size of the points in pixels (default: 2)
#' @param pt.shape Shape to use for the points (default: circle)
#' @param opacity Transparency level to use for the points on a 0-1 scale (default: 1)
#' @param palette.use Color palette to use.  Must be a palette available in the Paletteer package.  (default: 'Set1')
#' @param plot.height Plot height in pixels (default: 900)
#' @param plot.width Plot width in pixels (default: 900)
#' @param pt.info Meta.data columns to add to the hoverinfo popup. (default: ident)
#' @param legend Display legend? (default: TRUE)
#' @param legend.font.size Legend font size (default: 12)
#' @param plot.title Plot title (default: reduction.use)
#' @param plot.axes Display the major x, y, and z axes? (default: FALSE)
#' @param plot.grid Display the major unit tick marks? (default: FALSE)
#'
#' @import dplyr
#' @importFrom magrittr "%>%"
#' @importFrom plotly plot_ly
#' @importFrom plotly layout
#'
#' @return plotly object
#' @export
#'
#' @examples DimPlotly3D(seuratObj, group.by = "res.0.6", do.label = TRUE, show.arrow = FALSE)
#'
DimPlotly3D <- function(seuratObj,
                        group.by = "ident",
                        do.label = FALSE,
                        label.size = 12,
                        label.color = '000000',
                        show.arrow = FALSE,
                        do.return = FALSE,
                        pt.size = 2,
                        pt.shape = "circle",
                        opacity = 1,
                        reduction.use = "tsne",
                        dim.1 = 1,
                        dim.2 = 2,
                        dim.3 = 3,
                        palette.use = "Set1",
                        plot.height = 900,
                        plot.width = 900,
                        plot.title = NULL,
                        pt.info = NULL,
                        legend = TRUE,
                        legend.font.size = 12,
                        plot.grid = FALSE,
                        plot.axes = FALSE){

  df <- PrepDf(seuratObj,
               reduction.use,
               dim.1 = dim.1,
               dim.2 = dim.2,
               dim.3 = dim.3,
               ident = group.by)

  df <- PrepInfo(seuratObj = seuratObj,
                 pt.info = pt.info,
                 df = df)

  if (do.label) {
    df %>%
      dplyr::group_by(ident) %>%
      centers <- summarize(
        x = median(x = x),
        y = median(x = y),
        z = median(x = z)
      )
      labels <- list(x = centers$x,
                     y = centers$y,
                     z = centers$z,
                     text = centers$ident,
                     font = list(size = label.size)
      )
      compiled.labels = list()

      if (isTRUE(show.arrow)) {
        border.color = label.color
        bg.color = 'FFFFFF'
      } else {
        border.color = 'FFFFFF'
        bg.color = 'FFFFFF'
      }
      for (k in 1:length(unique(ident))) {
        tmp <- list(
          showarrow = show.arrow,
          x = labels$x[k],
          y = labels$y[k],
          z = labels$z[k],
          text = labels$text[k],
          font = list(size = label.size),
          bordercolor=border.color,
          bgcolor=bg.color,
          opacity=0.8)
        compiled.labels <- c(compiled.labels, list(tmp))
      }
  } else {
    compiled.labels <- NULL
  }

  if (is.null(plot.title)){
    plot.title <- reduction.use
  }

  pal <- PrepPalette(df = df,
                     palette.use = palette.use)

  p <- plot_ly(df,
               x = ~x,
               y = ~y,
               z = ~z,
               color = ~ident,
               colors = pal,
               marker = list(
                 symbol = pt.shape,
                 size = pt.size,
                 opacity = opacity,
                 mode = "markers"
               ),
               width = plot.width,
               height = plot.height,
               type = "scatter3d",
               mode = "markers",
               showlegend = legend
  ) %>%
    layout(
      title = plot.title,
      scene = list(
        aspectratio = list(x = 0,
                           y = 0,
                           z = -1),
        camera = list(
          center = list(x = 0,
                        y = 0,
                        z = 0
          ),
          eye = list(x = 2,
                     y = -1,
                     z = 0.5
          ),
          up = list(x = 1,
                    y = 0,
                    z = 0
          )
        ),
        dragmode = "turnable",
        xaxis = list(title = dim.axes[as.numeric(dim.1)],
                     type = "double",
                     showgrid = plot.grid,
                     visible = plot.axes
        ),
        yaxis = list(title = dim.axes[as.numeric(dim.2)],
                     type = "double",
                     showgrid = plot.grid,
                     visible = plot.axes
        ),
        zaxis = list(title = dim.axes[as.numeric(dim.3)],
                     type = "double",
                     showgrid = plot.grid,
                     visible = plot.axes
        ),
        annotations = compiled.labels
      )
    )

  if(!is.null(pt.info)){
    p <- p %>% add_markers(hoverinfo = "text",
                           hovertext = ~meta.info,
                           showlegend = FALSE
    )
  }

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
