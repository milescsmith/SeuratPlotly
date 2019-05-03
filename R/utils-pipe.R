#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL
#' Compound Assignment Pipe operator
#'
#' See \code{magrittr::\link[magrittr]{\%<>\%}} for details.
#'
#' @name %<>%
#' @rdname compound
#' @keywords internal
#' @export
#' @importFrom magrittr %<>%
#' @usage lhs \%<>\% rhs
NULL


#' @title NotIn
#' @name %nin%
#' @keywords internal
#' @export
#' @description Negative in
#' @usage lhs \%nin\% rhs
#'
#' @importFrom purrr compose
#'
`%nin%`` <- compose(`!`,`%in%`)
