library(endoU)
library(tidyverse)
library(fs)
library(usethis)

basedir <- "~/Dropbox (Hesselberth Lab)/Rachel_data/EndoU_project"

U6_tabs <- fs::path_join(c(basedir, "U6_bg"))
norm_tabs <- fs::path_join(c(basedir, "all_RNA_bg/exp_1"))

filenames <- list_names(norm_tabs)

library_sizes <- coverage_table(filenames) %>%
  dplyr::select(count, cell, virus, time) %>%
  dplyr::group_by(cell, virus, time) %>%
  dplyr::mutate(total_umi_reads = sum(count)) %>%
  dplyr::select(-count) %>%
  unique()

U6_cov_tbl <- coverage_table(list_names(U6_tabs))
U6_cov_tbl <- no_dinuc_complete(U6_cov_tbl)
U6_cov_tbl <- normalize_table(library_sizes, U6_cov_tbl)

usethis::use_data(U6_cov_tbl, compress = "xz", overwrite = TRUE)
