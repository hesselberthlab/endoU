#' Plot coverage
#'
#' @param data table object with count and position data
#' @param x variable for x-axis
#' @param y variable for y-axis
#' @param cell cell
#' @param virus virus
#' @param facet set to FALSE to generate graphs for a specifc cell and virus combination
#'
#'
#' @examples
#'
#' plot_coverage(viral_cov_small, pos, norm_count)
#'
#' plot_coverage(viral_cov_small, pos, norm_count, facet = FALSE, cell = "B6", virus = "ns2")
#'
#'
#' @export

plot_coverage <- function(data, x, y, cell = NULL, virus = NULL, facet = TRUE) {
  x <- enquo(x)
  y <- enquo(y)

  if (!is.null(cell)) {
    data <- dplyr::filter(data, cell == !!cell)
  }

  if (!is.null(virus)) {
    data <- dplyr::filter(data, virus == !!virus)
  }

  data <- data %>% dplyr::mutate(time = fct_relevel(time, "9", "12"))

  p <- ggplot(data, aes(x = (!!x/1000), y = !!y, fill = time, width = 0.025)) +
    geom_bar(stat = "identity", position = "identity") +
    scale_fill_endoU("time") +
    theme_cowplot() +
    ggtitle(paste(cell, virus, sep=" ")) +
    scale_x_continuous(breaks=seq(0, 30, 5)) +
    labs(
      x = "Position (kb)",
      y = "Normalized counts"
    )

  if (facet) {
    p <- p + facet_grid(virus ~ cell)
  }

  p
}


#' Plot fold change
#'
#' @param data table with log2 fold change data by cell or by virus
#' @param x variable for x-axis, based on a cell or virus fold change analysis
#' @param sample_comp variable to filter for speciifc fold change comparisions
#' @param facet set to FALSE to generate graphs for specifc fold change comparision
#'
#'
#' @examples
#'
#' plot_top_foldChange(data, cell)
#' plot_top_foldChange(data, cell, "MHVS_log2FC", facet = FALSE)
#'
#'
#' @export

plot_top_foldChange <- function(data, x, sample_comp = NULL, facet = TRUE) {
  x <- enquo(x)

  if (!is.null(sample_comp)) {
    data <- filter(data, sample_comp == !!sample_comp)
  }

  p <- ggplot(data, aes(x = (!!x), y = (-foldChange), color=factor(!!x))) +
    geom_quasirandom() +
    scale_color_endoU("general") +
    theme_cowplot() +
    ggtitle(paste(sample_comp)) +
    labs(
      y = "log2 fold change (WT/mutant)"
    )

  if (facet) {
    p <- p + facet_grid(~sample_comp)
  }

  p
}

#' Plot dinucleotides
#'
#' @param data table with percent dinucleotide data for cell, virus, and time
#' @param facet set to FALSE to generate graphs for specifc comparisions
#'
#'
#' @examples
#'
#' plot_dinuc(data)
#' plot_dinu(data, "B6", "MHVV, facet = FALSE)
#'
#'
#' @export
#'
#'

plot_dinuc <- function(data, x, y, cell = NULL, virus = NULL, facet = TRUE) {
  x <- enquo(x)
  y <- enquo(y)

  if (!is.null(cell)) {
    data <- filter(data, cell == !!cell)
  }

  if (!is.null(virus)) {
    data <- filter(data, virus == !!virus)
  }

  data <- data %>% mutate(time = fct_relevel(time, "9", "12"))

  p <- ggplot(data, aes(x = !!x, y = !!y, fill = time)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_endoU("time") +
    theme_cowplot() +
    ggtitle(paste(cell, virus, sep=" ")) +
    scale_x_discrete(limits=c("UA","UG", "UC", "UU","GA", "GG", "GC", "GU", "CA", "CG", "CC", "CU", "AA", "AG", "AC", "AU")) +
    labs(
      x = "Dinucleotide",
      y = "Percent total umi-reads"
    )

  if (facet) {
    p <- p + facet_grid(virus ~ cell)
  }

  p
}

#' Plot RNA capture by cell
#'
#' @param data table with column data containing percent or total reads by RNA type
#' @param y column name for RNA type of interest
#'
#'
#' @examples
#'
#' plot_RNA_capture(table, RNA_col)
#'
#' @export
#'
#'

plot_RNA_capture <- function(data, y) {

    y <- enquo(y)

    data <- data %>% mutate(time = fct_relevel(time, "9", "12"))

    p <- ggplot(data, aes(x = virus , y = !!y, fill = time)) +
    scale_fill_brewer(palette = 'Set1') +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_endoU("time") +
    theme_cowplot() +
    facet_grid( ~ cell) +
    ggtitle(y) +
    labs(
      y = "Percent total umi-reads"
    )
}

#' Plot RNA distribution by cell type
#'
#' @param data table with column data containing percent reads by RNA type
#' @param time time hpi to filter on
#' @param cols_plot column numbers identifying percent RNA captured for RNAs of interest
#' @param cell set to null, can provide cell type of interest
#' @param facet set to FALSE if providing specific cell type
#'
#'
#' @examples
#'
#' plot_RNAbyCell(final_table, "12", 13:20, "B6", facet = FALSE)
#'
#' plot_RNAbyCell(final_table, "12", 13:20)
#'
#' @export
#'

plot_RNAbyCell <- function(data, time, cols_plot, cell = NULL, facet = TRUE) {

  if (!is.null(cell)) {
    data <- filter(data, cell == !!cell)
  }

    data <- filter(data, time == !!time)

    res <- data %>%
      gather(RNA_type, percent_total_reads, cols_plot) %>%
      separate(RNA_type, into = c('RNA_type'), sep = '_')

    p <- ggplot(res, aes(x = RNA_type , y = percent_total_reads)) +
      geom_bar(aes(fill = virus), stat = "identity", position = 'dodge') +
      scale_fill_endoU("general") +
      theme_cowplot() +
      ggtitle(paste(cell, paste(time, "hpi", sep = "_"), sep=" ")) +
      labs(
        y = "Percent total umi-reads"
      )

    if (facet) {
      p <- p + facet_grid(~ cell)
    }

    p
}

#Palettes ----------------------------------------------------------------------
#' Function to extract endoU colors as hex codes
#'
#' @examples
#'
#' geom_point(color = endoU_cols("sky blue))

endoU_colors <- c(
  `black`        = "#000000",
  `orange`      = "#E69F00",
  `sky blue`       = "#56B4E9",
  `bluish green`     = "#009E73",
  `yellow`     = "#F0E442",
  `blue` = "#0072B2",
  `vermillion`  = "#D55E00",
  `reddish purple` = "#CC79A7",
  `dark grey`  = "#8c8c8c")


endoU_cols <- function(...) {
  cols <- c(...)

  if (is.null(cols))
    return (endoU_colors)

  endoU_colors[cols]
}

#' Return function to interpolate a endoU color palette
#'
#' @param palette Character name of palette in endoU_palettes
#' @param reverse Boolean indicating whether the palette should be reversed
#' @param ... Additional arguments to pass to colorRampPalette()
#'

endoU_palettes <- list(
  `time`  = endoU_cols("orange", "sky blue"),

  `cool`  = endoU_cols("sky blue", "bluish green", "blue"),

  `hot`   = endoU_cols("yellow", "orange", "vermillion, reddish purple"),

  `all` = endoU_cols("black", "orange", "sky blue", "bluish green", "yellow", "blue", "vermillion", "reddish purple", "dark grey"),

  `general` = endoU_cols("bluish green", "yellow", "blue", "vermillion", "black", "reddish purple", "dark grey"),

  `grey`  = endoU_cols("black", "dark grey")
)


endoU_pal <- function(palette = "general", reverse = FALSE, ...) {
  pal <- endoU_palettes[[palette]]

  if (reverse) pal <- rev(pal)

  colorRampPalette(pal, ...)
}

#' Color scale constructor for endoU colors
#'
#' @param palette Character name of palette in endoU_palettes
#' @param discrete Boolean indicating whether color aesthetic is discrete or not
#' @param reverse Boolean indicating whether the palette should be reversed
#' @param ... Additional arguments passed to discrete_scale() or
#'            scale_color_gradientn(), used respectively when discrete is TRUE or FALSE
#'
#' @examples scale_color_endoU("grey")
#'
#' @export
#'

scale_color_endoU <- function(palette = "general", discrete = TRUE, reverse = FALSE, ...) {
  pal <- endoU_pal(palette = palette, reverse = reverse)

  if (discrete) {
    discrete_scale("color", paste0("endoU_", palette), palette = pal, ...)
  } else {
    scale_color_gradientn(colours = pal(256), ...)
  }
}

#' Fill scale constructor for endoU colors
#'
#' @param palette Character name of palette in endoU_palettes
#' @param discrete Boolean indicating whether color aesthetic is discrete or not
#' @param reverse Boolean indicating whether the palette should be reversed
#' @param ... Additional arguments passed to discrete_scale() or
#'            scale_fill_gradientn(), used respectively when discrete is TRUE or FALSE
#'
#' @examples scale_fill_endoU("grey")
#'
#' @export
#'

scale_fill_endoU <- function(palette = "general", discrete = TRUE, reverse = FALSE, ...) {
  pal <- endoU_pal(palette = palette, reverse = reverse)

  if (discrete) {
    discrete_scale("fill", paste0("endoU_", palette), palette = pal, ...)
  } else {
    scale_fill_gradientn(colours = pal(256), ...)
  }
}

