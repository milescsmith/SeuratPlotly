#' GOBubblePlotly
#'
#' Produces a Bubble Plot for the genes of a given GO term.
#'
#' @param object Seurat object
#' @param go_term Gene Ontology term identifier (i.e. GO:0046774)
#' @param grouping_var Factor by which to group cells.  (default: ident)
#' @param plot_height Plot height in pixels. (default: 900)
#' @param gene_filter
#' @param ...
#' @param plot_width Plot width in pixels. (default: 900)
#'
#' @importFrom dplyr select distinct filter
#'
#' @return
#' @export
#'
#' @examples
GOBubblePlotly <- function(object,
                           go_term,
                           grouping_var = "ident",
                           plot_height = 900,
                           plot_width = 900,
                           gene_filter = NULL,
                           ...){

  if(missing(go_term)){
    stop("No GO term supplied")
  }
  go_genes_to_plot <- retrieveGO(go_term) %>%
                      select(hgnc_symbol) %>%
                      distinct() %>%
  if (!is.null(gene_filter)){
    filter(hgnc_symbol %in% gene_filter)
  }

  if(length(go_genes_to_plot) > 0){
    BubblePlotly(object,
                 features_plot = go_genes_to_plot$hgnc_symbol,
                 grouping_var = grouping_var,
                 plot_height = plot_height,
                 plot_width = plot_width,
                 ...)
  } else {
    print("No genes for that term are expressed in the dataset.")
  }

}
