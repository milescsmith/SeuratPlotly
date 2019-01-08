#' GOBubblePlotly
#'
#' Produces a Bubble Plot for the genes of a given GO term.
#'
#' @param seuratObj Seurat object
#' @param go_term Gene Ontology term identifier (i.e. GO:0046774)
#' @param group.by Factor by which to group cells.  (default: ident)
#' @param plot.height Plot height in pixels. (default: 900)
#' @param plot.width Plot width in pixels. (default: 900)
#' @param filter A list of gene names to filter the GO term members against. (default: seuratObj@var.genes)
#' @param ...options to pass to BubblePlotly
#'
#' @import dplyr
#' @importFrom magrittr "%>%"
#'
#' @return
#' @export
#'
#' @examples
GOBubblePlotly <- function(seuratObj,
                           go_term,
                           group.by = "ident",
                           plot.height = 900,
                           plot.width = 900,
                           filter = seuratObj@var.genes,
                           ...){

  go_genes_to_plot <- retrieveGO(go_term) %>%
                      select(hgnc_symbol) %>%
                      distinct() %>%
                      filter(hgnc_symbol %in% filter)

  if(length(go_genes_to_plot) > 0){
    BubblePlotly(seuratObj,
                 genes.plot = go_genes_to_plot$hgnc_symbol,
                 group.by = group.by,
                 plot.height = plot.height,
                 plot.width = plot.width,
                 ...)
  } else {
    print("No genes for that term are expressed in the dataset.")
  }

}
