#' GOBubblePlotly
#'
#' Produces a Bubble Plot for the genes of a given GO term.
#'
#' @param object Seurat object
#' @param go_term Gene Ontology term identifier (i.e. GO:0046774)
#' @param grouping Factor by which to group cells.  (default: ident)
#' @param plot_height Plot height in pixels. (default: 900)
#' @param plot_width Plot width in pixels. (default: 900)
#' @param filter A list of gene names to filter the GO term members against. (default: object@var.genes)
#' @param ...options to pass to BubblePlotly
#'
#' @import dplyr
#' @importFrom magrittr "%>%"
#'
#' @return
#' @export
#'
#' @examples
GOBubblePlotly <- function(object,
                           go_term,
                           grouping = "ident",
                           plot_height = 900,
                           plot_width = 900,
                           filter = object@var.genes,
                           ...){

  go_genes_to_plot <- retrieveGO(go_term) %>%
                      select(hgnc_symbol) %>%
                      distinct() %>%
                      filter(hgnc_symbol %in% filter)

  if(length(go_genes_to_plot) > 0){
    BubblePlotly(object,
                 genes.plot = go_genes_to_plot$hgnc_symbol,
                 grouping = grouping,
                 plot_height = plot_height,
                 plot_width = plot_width,
                 ...)
  } else {
    print("No genes for that term are expressed in the dataset.")
  }

}
