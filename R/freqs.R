#' Load frequencies from a file
#'
#' @param x tsv file with frequencies
#' @examples
#' # path to file
#' load_freqs(x)
#'
#' @export
load_freqs <- function(x) {
  x <- readr::read_tsv(
    col_names = c("dinuc", "count")
  )

  x <- mutate(freq = count / sum(count))
  x <- rename(
    x,
    mhv_count = count,
    mhv_freq = freq
  )

  x
}

#' @export
calc_freqs <- function(x) {
  df <- readr::read_tsv(
    x,
    col_names = c("chrom", "start", "end", "count", "normalized_count", "dinuc")
  )

  df <- left_join(df)
}
