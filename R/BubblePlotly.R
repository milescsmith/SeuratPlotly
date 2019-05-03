# TODO: add more ways to sort/slice/dice!
#' BubblePlotly
#'
#' Plot average expression levels and proportion expressing a feature arranged by given grouping_var_var.
#' Color intensity of the markers indicates expression levels and size of the marker indicates
#' the proportion of the group expressing the gene above a given threshold. By default, genes are arranged
#' along the x-axis and groups along the y-axis.
#'
#' @param object Seurat object
#' @param genes_plot A list of genes to plot for each group.
#' @param colors_use Color palette to use.  Palettes from RColorBrewer and viridis. Default: Reds
#' @param dot_min Minimium marker size, in pixels. Default: 0
#' @param dot_scale Factor by which to scale markers. Default: 2
#' @param grouping_var Factor by which to group cells.  Default: ident
#' @param legend Display legend. (currently nonfunctional) Default TRUE
#' @param return Return the plot object instead of displaying it Default: FALSE
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
#' @param flip Swap the x- and y-axes so that genes are on the y-axis and groups along the x-axis. Default: false
#' @param alphabetize Alphabetize the display order of genes. Default: TRUE
#' @param pct_expr_thres Hide a gene if there is no group that expresses at or above this percentage. Default: 0
#' @param return Return the generated data frame underlying the graph. Default: FALSE
#' @param use_scaled Use scaled data. Default: FALSE
#' @param use_raw Use raw data. Default: FALSE
#' @param show_zeros Display comparison groups for which there there was no data? Default: FALSE
#' @param add_group Add a null comparison group.  Will be initialized with zeros.
#' @param y_label_order List of labels to use for the comparison groups.  Will be displayed in order of their index.
#'
#' @importFrom dplyr filter summarise group_by mutate inner_join ungroup distinct arrange mutate_if
#' @importFrom tidyr gather
#' @importFrom tibble rownames_to_column add_row
#' @importFrom compositions normalize
#' @importFrom plotly plot_ly layout
#' @importFrom gtools mixedsort
#'
#' @return
#' @export
#'
#' @examples
BubblePlotly <- function (object,
                          genes_plot,
                          assay_use = "RNA",
                          slot_use = "data",
                          grouping_var = "ident",
                          colors_use = "Blues",
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

  #screen out any genes that are not in our dataset and print them
  original_genes_to_plot <- genes_plot
  genes_plot <- genes_plot[which(genes_plot %in% rownames(object))]
  not_found <- original_genes_to_plot[original_genes_to_plot %nin% genes_plot]
  message(not_found)

  #let the user alphabetize the list of genes.  makes it easier to find a particular gene in a big table
  #TODO: add more sort options
  if (isTRUE(alphabetize)){
    genes_plot <- mixedsort(genes_plot, decreasing = TRUE)
  }

  df <- GetFeatureValues(object = object,
                         features = genes_plot,
                         assay_use = assay_use,
                         slot_use = slot_use,
                         bins = bins)

  # Add 0 data for each gene that was not detected or failed to pass QC
  if (isTRUE(show_zeros)){
    for(i in not_found){
      df[,i] <- 0
    }
  }

  ident <- GetFeatureValues(object, features = grouping_var)
  colnames(ident)[1] <- "ident"
  df <- inner_join(df, rownames_to_column(ident, "cell"))

  df <- df %>% gather(key = genes_plot,
                      value = expression,
                      -c(cell, ident))

  df <- df %>% group_by(ident, genes_plot) %>%
    summarise(avg.exp = mean(x = expression),
              pct.exp = PercentAbove(x = expression, threshold = 0),
              n = n())

  if(!is.null(pct_expr_thres)){
    df <- df %>% group_by(genes_plot) %>% filter(max(pct.exp) > pct_expr_thres)
  }

  df <- df %>% ungroup() %>% group_by(genes_plot) %>%
    mutate(avg.exp.scale = replace_na(normalize(avg.exp), 0))

  if (isTRUE(show_zeros)){
    df[["genes_plot"]] <- factor(x = df[["genes_plot"]],
                                      levels = rev(x = sub(pattern = "-",
                                                           replacement = ".",
                                                           x = c(genes_plot, not_found))))
  } else {
    df[["genes_plot"]] <- factor(x = df[["genes_plot"]],
                                 levels = rev(x = sub(pattern = "-",
                                                      replacement = ".",
                                                      x = genes_plot)))
  }

  df[["pct.exp"]][df[["pct.exp"]] < dot_min] <- NA

  if(!is.null(add_group)){
    for(i in add_group){
      for(j in genes_plot){
        df <- df %>%
          do(add_row(.,ident = i,
                     genes_plot = j,
                     avg.exp = 0,
                     pct.exp = 0,
                     n = 0,
                     avg.exp.scale=0)) %>%
          distinct()
      }
    }
  }

  pal <- PrepQualitativePalette(bins = bins, palette_use = colors_use)

  df <- df %>% arrange(ident) %>% ungroup() %>% mutate_if(is.factor, list(~as.character(.)))

  genes_plot <- df$genes_plot

  if(flip){
    ax = ~ident
    ay = ~genes_plot
  } else {
    ax = ~genes_plot
    ay = ~ident
  }

  p <- plot_ly(data = df,
               x = ~ident,
               y = ~genes_plot,
               # x = ax,
               # y = ay,
               hoverinfo = 'text',
               text = ~paste('Group:', ident,
                             '<br>Gene:', genes_plot,
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
           xaxis = list(title = "Gene", tickangle = x_label_rotation, tickfont = list(size = x_font_size)),
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
  # print(length(unique(df$id)))
  # print(length(unique(df$genes_plot)))
}
