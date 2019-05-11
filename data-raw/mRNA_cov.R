library(endoU)
library(tidyverse)
library(fs)
library(usethis)

basedir <- "~/Dropbox (Hesselberth Lab)/Rachel_data/EndoU_project"

mRNA_tabs <- fs::path_join(c(basedir, "new_mRNA"))
norm_tabs <- fs::path_join(c(basedir, "all_RNA_bg/exp_1"))

filenames <- list_names(norm_tabs)

library_sizes <- coverage_table(filenames) %>%
  dplyr::select(count, cell, virus, time) %>%
  dplyr::group_by(cell, virus, time) %>%
  dplyr::mutate(total_umi_reads = sum(count)) %>%
  dplyr::select(-count) %>%
  unique()

mRNA_cov_tbl <- coverage_table(list_names(mRNA_tabs))
mRNA_cov_tbl <- chrom_normalize_table(library_sizes, mRNA_cov_tbl)

usethis::use_data(mRNA_cov_tbl, compress = "xz", overwrite = TRUE)
