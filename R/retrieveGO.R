#' retrieveGO
#'
#' Retrieve the HUGO Gene Nomenclature Committee names associated with a GeneOntology term.
#'
#' @param term GO ID to retreive names for.
#' @param mart Biomart database to use. Default: 'ensembl'
#' @param dataset Biomart dataset to use. Default: 'hsapiens_gene_ensembl'
#'
#' @importFrom biomaRt useMart getBM
#'
#' @return A list of HCGN names.
#' @export
#'
#' @examples genes_to_plot <- retrieveGO('GO:0046774'
#'
retrieveGO <- function(term, mart = 'ensembl', dataset='hsapiens_gene_ensembl'){
  gene_data <- getBM(attributes=c('hgnc_symbol', 'ensembl_transcript_id', 'go_id'),
                     filters = 'go',
                     values = term,
                     mart = useMart(mart,dataset=dataset)
                     )
  return(gene_data)
}
