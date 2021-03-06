library(endoU)
library(tidyverse)
library(fs)
library(usethis)

basedir <- "~/Dropbox (Hesselberth Lab)/Rachel_data/EndoU_project"

r28S_tabs <- fs::path_join(c(basedir, "28S/exp1/new"))
norm_tabs <- fs::path_join(c(basedir, "all_RNA_bg/exp_1"))
dinuc_counts <- fs::path_join(c(basedir, "dinuc_byposition/r28s.dinuc.txt"))
dinuc_counts1 <- fs::path_join(c(basedir, "dinuc_byposition/r28s.dinuc_+1.txt"))

filenames <- list_names(norm_tabs)

library_sizes <- coverage_table(filenames) %>%
  dplyr::select(count, cell, virus, time) %>%
  dplyr::group_by(cell, virus, time) %>%
  dplyr::mutate(total_umi_reads = sum(count)) %>%
  dplyr::select(-count) %>%
  unique()

r28S_cov_tbl <- coverage_table(list_names(r28S_tabs))
r28S_cov_tbl <- complete_table(r28S_cov_tbl, dinuc_counts)
r28S_cov_tbl <- normalize_table(library_sizes, r28S_cov_tbl)

r28S_cov_tbl1 <- coverage_table(list_names(r28S_tabs))
r28S_cov_tbl1 <- complete_table(r28S_cov_tbl1, dinuc_counts1)
r28S_cov_tbl1 <- normalize_table(library_sizes, r28S_cov_tbl1)

usethis::use_data(r28S_cov_tbl, compress = "xz", overwrite = TRUE)
usethis::use_data(r28S_cov_tbl1, compress = "xz", overwrite = TRUE)
