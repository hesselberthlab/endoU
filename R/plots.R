#' Plot coverage
#'
#' @param data table object with count and position data
#' @param x variable for x-axis
#' @param y variable for y-axis
#' @param cell cell
#' @param virus virus
#'
#' @examples
#' plot_coverage(
#'   viral_cov_small,
#'   start, norm_count,
#'   cell = "B6", virus = "ns2"
#' )
#'
#' @export
plot_coverage <- function(data, x, y, cell = NULL, virus = NULL) {
  x <- enquo(x)
  y <- enquo(y)

  if (!is.null(cell)) {
    data <- filter(data, cell == !!cell)
  }

  if (!is.null(virus)) {
    data <- filter(data, virus == !!virus)
  }

  ggplot(data, aes(x = !!x, y = !!y, fill = time, width = 25)) +
    geom_bar(stat = "identity", position = "identity") +
    scale_fill_brewer(palette = "Set1") +
    theme_cowplot() +
    labs(
      x = "Position (kb)",
      y = "Normalized counts"
    )

}



