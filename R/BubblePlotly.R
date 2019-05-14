#' @title BubblePlotly
#'
#' @description Plot average expression levels and proportion expressing a feature arranged by given grouping_var.
#' Color intensity of the markers indicates expression levels and size of the marker indicates
#' the proportion of the group expressing the feature above a given threshold. By default, features are arranged
#' along the x-axis and groups along the y-axis.
#'
#' @param object Seurat object
#' @param features_plot Features to display values for. If trying to retrieve values from an assay
#' that is not the current DefaultAssay, prefix the feature name with 'assay_'.
#' Works with anything \code{Seurat::\link[Seurat]{FetchData}} can retrieve.
#' @param assay Assay to pull data from.  Defaults to the object's DefaultAssay
#' @param slot Slot to pull data from.  Default: 'data'
#' @param grouping_var Factor by which to group cells.  Default: ident
#' @param colors Color palette to use. Default: Reds
#' @param dot_min Minimium marker size, in pixels. Default: 0
#' @param dot_scale Factor by which to scale markers. Default: 2
#' @param plot_legend Display legend? Default: FALSE
#' @param x_label_rotation Angle by which to rotate the x-axis labels, in degrees relative to horizontally aligned text. Default -45°
#' @param plot_height Plot height in pixels. Default: 900
#' @param plot_width Plot width in pixels. Default: 900
#' @param x_font_size Size of the x-axis titles. Default: 10
#' @param y_font_size Size of the y-axis titles. Default: 10
#' @param title_font_size Size of the plot title. Default: 12
#' @param legend_text_size Size of the legend text. Default: 10
#' @param opacity Transparency level to use for the points, on a 0-1 scale Default: 1
#' @param plot_title  Display title with the name of the feature? Default TRUE
#' @param bins Number of bins to use in dividing expression levels. Default: 10
#' @param flip Swap the x- and y-axes so that features are on the y-axis and groups along the x-axis. Default: false
#' @param alphabetize Alphabetize the display order of features Default: TRUE
#' @param pct_expr_thres Hide a feature if there is no group that expresses at or above this percentage. Default: 0
#' @param return Return the generated data frame underlying the graph. Default: FALSE
#' @param show_zeros Display comparison groups for which there there was no data? Default: FALSE
#' @param add_group Add a null comparison group.  Will be initialized with zeros.
#' @param y_label_order List of labels to use for the comparison groups.  Will be displayed in order of their index.
#'
#' @importFrom dplyr filter summarise group_by mutate inner_join
#' ungroup distinct arrange mutate_if do
#' @importFrom tidyr gather replace_na
#' @importFrom tibble rownames_to_column add_row
#' @importFrom compositions normalize
#' @importFrom plotly plot_ly layout
#' @importFrom gtools mixedsort
#'
#' @return
#' @export
#'
BubblePlotly <- function (object,
                          features_plot,
                          assay = NULL,
                          slot = "data",
                          grouping_var = "ident",
                          colors = "Blues",
                          dot_min = 0,
                          dot_scale = 2,
                          plot_legend = FALSE,
                          return = FALSE,
                          x_label_rotation = -45,
                          plot_width = 600,
                          plot_height = 600,
                          x_font_size = 10,
                          y_font_size = 10,
                          title_font_size = 12,
                          legend_text_size = 10,
                          opacity = 1,
                          plot_title = NULL,
                          bins = 50,
                          flip = FALSE,
                          alphabetize = TRUE,
                          pct_expr_thres = NULL,
                          show_zeros = FALSE,
                          add_group = NULL,
                          y_label_order = NULL){

  #globalbindings stuff
  cell <- NULL
  pct.exp <- NULL
  avg.exp <- NULL
  n <- NULL

  #screen out any features that are not in our dataset and print them
  original_features_to_plot <- features_plot
  features_plot <- features_plot[which(features_plot %in% rownames(object))]
  not_found <- original_features_to_plot[original_features_to_plot %nin% features_plot]
  message(not_found)

  #let the user alphabetize the list of features.  makes it easier to find a particular feature in a big table
  #TODO: add more sort options
  if (isTRUE(alphabetize)){
    features_plot <- mixedsort(features_plot, decreasing = TRUE)
  }

  df <- GetFeatureValues(object = object,
                         features = features_plot,
                         assay = assay,
                         slot = slot,
                         bins = bins)

  # Add 0 data for each feature that was not detected or failed to pass QC
  if (isTRUE(show_zeros)){
    for(i in not_found){
      df[,i] <- 0
    }
  }

  ident <- GetFeatureValues(object, features = grouping_var)
  colnames(ident)[1] <- "ident"
  df <- inner_join(df, rownames_to_column(ident, "cell"))

  df <- df %>% gather(key = features_plot,
                      value = expression,
                      -c(cell, ident))

  df <- df %>% group_by(ident, features_plot) %>%
    summarise(avg.exp = mean(x = expression),
              pct.exp = PercentAbove(x = expression, threshold = 0),
              n = n())

  if(!is.null(pct_expr_thres)){
    df <- df %>% group_by(features_plot) %>% filter(max(pct.exp) > pct_expr_thres)
  }

  df <- df %>% ungroup() %>% group_by(features_plot) %>%
    mutate(avg.exp.scale = replace_na(normalize(avg.exp), 0))

  if (isTRUE(show_zeros)){
    df[["features_plot"]] <- factor(x = df[["features_plot"]],
                                      levels = rev(x = sub(pattern = "-",
                                                           replacement = ".",
                                                           x = c(features_plot, not_found))))
  } else {
    df[["features_plot"]] <- factor(x = df[["features_plot"]],
                                 levels = rev(x = sub(pattern = "-",
                                                      replacement = ".",
                                                      x = features_plot)))
  }

  df[["pct.exp"]][df[["pct.exp"]] < dot_min] <- NA

  if(!is.null(add_group)){
    for(i in add_group){
      for(j in features_plot){
        df <- df %>%
          do(add_row(.,ident = i,
                     features_plot = j,
                     avg.exp = 0,
                     pct.exp = 0,
                     n = 0,
                     avg.exp.scale=0)) %>%
          distinct()
      }
    }
  }

  pal <- PrepQuantitativePalette(bins = bins, palette = colors)

  df <- df %>% arrange(ident) %>% ungroup() %>% mutate_if(is.factor, list(~as.character(.)))

  features_plot <- df$features_plot

  if(flip){
    ax = ~ident
    ay = ~features_plot
  } else {
    ax = ~features_plot
    ay = ~ident
  }

  p <- plot_ly(data = df,
               x = ~ident,
               y = ~features_plot,
               # x = ax,
               # y = ay,
               hoverinfo = 'text',
               text = ~paste('Group:', ident,
                             '<br>feature:', features_plot,
                             '<br>Avg normalized expression:', round(avg.exp.scale,2),
                             '<br>% expressing:', round(pct.exp*100,2)),
               type = 'scatter',
               mode = 'markers',
               color = ~avg.exp.scale,
               colors = pal,
               cmin = pal[1],
               cmax = pal[100],
               size = ~pct.exp,
               sizes = c(0,10)*dot_scale,
               marker = list(opacity = opacity,
                             symbol = 'circle',
                             sizemode = 'diameter'),
               showlegend = FALSE,
               width = plot_width,
               height = plot_height) %>%
    layout(title = plot_title,
           titlefont = list(size = title_font_size),
           #font = list(size = font.size),
           xaxis = list(title = "feature", tickangle = x_label_rotation, tickfont = list(size = x_font_size)),
           yaxis = list(title = "Cluster", tickfont = list(size = y_font_size), type = "category", categoryorder = "array", categoryarray = y_label_order),
           margin = list(b = 100, l = 250)
    )

  if (return) {
    return(p)
  } else if(isTRUE(return))
  {
    return(df)
  } else {
    p
  }
}
