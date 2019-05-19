library(endoU)
library(tidyverse)
library(fs)
library(usethis)

basedir <- "~/Dropbox (Hesselberth Lab)/Rachel_data/EndoU_project"
norm_tabs <- fs::path_join(c(basedir, "all_RNA_bg/exp_1"))
viral_tabs <- fs::path_join(c(basedir, "viral_bg/pos"))
dinuc_positions <- fs::path_join(c(basedir, "dinuc_byposition/mhv.dinuc.sense.txt"))

filenames <- list_names(norm_tabs)

library_sizes <- coverage_table(filenames) %>%
  dplyr::select(count, cell, virus, time) %>%
  dplyr::group_by(cell, virus, time) %>%
  dplyr::mutate(total_umi_reads = sum(count)) %>%
  dplyr::select(-count) %>%
  unique()

antiviral_cov_tbl <- coverage_table(list_names(viral_tabs))
antiviral_cov_tbl <- complete_table(antiviral_cov_tbl, dinuc_positions)
antiviral_cov_tbl <- normalize_table(library_sizes, antiviral_cov_tbl)

usethis::use_data(antiviral_cov_tbl, compress = "xz", overwrite = TRUE)
