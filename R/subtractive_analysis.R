#' Normalize count data using the signal that occurs in the abscence of RNaseL
#'
#' @param x viral_cov_table
#'
#' @examples
#'
#' RNaseL_sub_tbl <- RNaseL_sub(viral_cov_table)
#'
#'
#' @export
#'

RNaseL_sub <- function(x) {
  dplyr::select(x, -count, -dinuc) %>%
  spread(cell, norm_count)
}


