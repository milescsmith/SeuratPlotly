#' PercentAbove
#'
#' Return the percentage of a list that is above a threshold
#'
#' @param x A list of numeric values.
#' @param threshold A numeric threshold value.
#'
#' @return double
#' @export
#'
#' @examples
PercentAbove <- function(x, threshold){
  return(length(x = x[x > threshold]) / length(x = x))
}
