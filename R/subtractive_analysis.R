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




#spread(cell, norm_count)
  subtractive_neg[is.na(subtractive_neg)] <- 0

subtractive_neg <- subtractive_neg %>%
  mutate(B6_RNaseLdep = B6 - RNaseL) %>%
  mutate(IFNAR_RNaseLdep = IFNAR - RNaseL) %>%
  dplyr::select(start:end, virus:time, B6_RNaseLdep:IFNAR_RNaseLdep) %>%
  gather(cell, normalized_count, B6_RNaseLdep:IFNAR_RNaseLdep)
