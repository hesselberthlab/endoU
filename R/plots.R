#' Exploratory plots
#'
#' @param data table object with count and position data
#' @param x variable for x-axis
#' @param y variable for y-axis
#'
#' @export
plot_depth_multi <- function(data, x, y) {
  x <- enquo(x)
  y <- enquo(y)

  ggplot(data, aes(x = !!x, y = !!y, fill = time, width = 25)) +
    geom_bar(stat = "identity", position = "identity") +
    scale_fill_brewer(palette = "Set1") +
    theme_cowplot() +
    facet_grid(virus ~ cell) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
    labs(x="Position", y="Normalized counts") +
    scale_x_continuous(breaks=seq(0, 30000, 5000))
}

#' Plot depth for a single variable set
#'
#' @inheritParams plot_depth_multi
#' @param cell cell
#' @param virus virus
#'
#' @examples
#' plot_depth_single(
#'   viral_cov_tbl,
#'   start, norm_count,
#'   B6, ns2
#' )
#'
#' @export
plot_depth_single <- function(data, x, y, cell, virus) {
  x <- enquo(x)
  y <- enquo(y)
  cell <- enquo(cell)
  virus <- enquo(virus)

  data <- filter(data, cell == !!cell & virus == !!virus) %>%
    dplyr::mutate(start == (start/1000))

  ggplot(data, aes(x = !!x, y = !!y, fill = time, width = 25)) +
    geom_bar(stat = "identity", position = "identity") +
    scale_fill_brewer(palette = "Set1") +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
    labs(x="Position (kb)", y="Normalized counts")
}

 #cells <- c("B6", "IFNAR", "RNASEL")
 #virus <- c("MHVS", "MHVV", "mock", "ns2", "nsp15")
 #combos <- tidyr::crossing(cells, virus)
