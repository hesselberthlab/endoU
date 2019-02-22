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
#'   pos, norm_count,
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

  data <- data %>% mutate(time = fct_relevel(time, "9", "12"))

  ggplot(data, aes(x = (!!x/1000), y = !!y, fill = time, width = 0.025)) +
    geom_bar(stat = "identity", position = "identity") +
    scale_fill_brewer(palette = "Set1") +
    theme_cowplot() +
    ggtitle(paste(cell, virus, sep="_")) +
    scale_y_continuous(limits = c(0.0, 0.2)) +
    scale_x_continuous(breaks=seq(0, 30, 5)) +
    labs(
      x = "Position (kb)",
      y = "Normalized counts"
    )

}
