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
    data <- filter(data, cell == !!cell)
  }

  if (!is.null(virus)) {
    data <- filter(data, virus == !!virus)
  }

  data <- data %>% mutate(time = fct_relevel(time, "9", "12"))

  p <- ggplot(data, aes(x = (!!x/1000), y = !!y, fill = time, width = 0.025)) +
    geom_bar(stat = "identity", position = "identity") +
    scale_fill_OkabeIto() +
    theme_cowplot() +
    ggtitle(paste(cell, virus, sep=" ")) +
    scale_y_continuous(limits = c(0.0, 0.21)) +
    scale_x_continuous(breaks=seq(0, 30, 5)) +
    labs(
      x = "Position (kb)",
      y = "Normalized counts"
    ) +
    theme(
      axis.title.x = element_text(size=11),
      axis.title.y = element_text(size=11),
      plot.title = element_text(size=12, face="bold", hjust=0))

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
    scale_color_OkabeIto() +
    theme_cowplot() +
    ggtitle(paste(sample_comp)) +
    #geom_hline(yintercept = 0) +
    labs(
      y = "log2 fold change (WT/endoU mutant)"
    ) +
    theme(
      axis.title.x = element_text(size=11),
      axis.title.y = element_text(size=11),
      plot.title = element_text(size=12, face="bold", hjust=0))

  if (facet) {
    p <- p + facet_grid(~sample_comp)
  }

  p
}

#' Plot dinucleotides
#'
#' @param data table with percent dinucleotide data for cell, virus, and time
#' @param facet set to FALSE to generate graphs for specifc fold change comparision
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
    scale_fill_OkabeIto() +
    theme_cowplot() +
    ggtitle(paste(cell, virus, sep=" ")) +
    scale_x_discrete(limits=c("UA","UG", "UC", "UU","GA", "GG", "GC", "GU", "CA", "CG", "CC", "CU", "AA", "AG", "AC", "AU")) +
    labs(
      x = "Dinucleotide",
      y = "Percent total umi-reads"
    ) +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1, size=11),
      axis.title.x = element_text(size=11),
      axis.title.y = element_text(size=11),
      plot.title = element_text(size=12, face="bold", hjust=0))

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

    p <- ggplot(data, aes(x = virus , y = !!y, fill = time)) +
    scale_fill_brewer(palette = 'Set1') +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_OkabeIto() +
    theme_cowplot() +
    facet_grid( ~ cell) +
    ggtitle(y) +
    labs(
      y = "Percent total umi-reads"
    ) +
      theme(
        axis.text.x = element_text(angle = 60, hjust = 1, size=11),
        axis.title.x = element_text(size=11),
        axis.title.y = element_text(size=11),
        plot.title = element_text(size=12, face="bold", hjust=0))
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
#' plot_RNAbyCell(final_table, "RNaseL", "12", 13:20)
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
      scale_fill_OkabeIto() +
      theme_cowplot() +
      ggtitle(paste(cell, paste(time, "hpi", sep = "_"), sep=" ")) +
      labs(
        y = "Percent total umi-reads"
      ) +
      theme(
        axis.text.x = element_text(angle = 60, hjust = 1, size=11),
        axis.title.x = element_text(size=11),
        axis.title.y = element_text(size=11),
        plot.title = element_text(size=12, face="bold", hjust=0))


    if (facet) {
      p <- p + facet_grid(~ cell)
    }

    p
}

