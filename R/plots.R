#' Exploratory plots
#'
#' @param x table object with count and position data
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

plot_depth_single <- function(data, x, y, cell, virus) {
  x <- enquo(x)
  y <- enquo(y)

  data <- filter(data, cell == cell & virus == virus) %>%
    dplyr::mutate(start == (start/1000))

  ggplot(data, aes(x = !!x, y = !!y, fill = time, width = 25)) +
    geom_bar(stat = "identity", position = "identity") +
    scale_fill_brewer(palette = "Set1") +
    theme_cowplot() +
    # facet_grid(virus ~ cell) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
    labs(x="Position (kb)", y="Normalized counts")
    #scale_x_continuous(breaks=seq(0, 30000, 5000))
}

 #cells <- c("B6", "IFNAR", "RNASEL")
 #virus <- c("MHVS", "MHVV", "mock", "ns2", "nsp15")
 #combos <- tidyr::crossing(cells, virus)
