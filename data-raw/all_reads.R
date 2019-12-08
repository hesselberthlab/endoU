library(endoU)
library(tidyverse)
library(fs)
library(usethis)

basedir <- "~/Dropbox (Hesselberth Lab)/Rachel_data/EndoU_project"
norm_tabs <- fs::path_join(c(basedir, "all_RNA_bg/exp_2"))

filenames <- list_names(norm_tabs)

library_sizes <- coverage_table(filenames) %>%
  dplyr::select(count, cell, virus, time) %>%
  dplyr::group_by(cell, virus, time) %>%
  dplyr::mutate(total_umi_reads = sum(count)) %>%
  dplyr::select(-count) %>%
  unique()

usethis::use_data(library_sizes, compress = "xz", overwrite = TRUE)

basedir <- "~/Dropbox (Hesselberth Lab)/Rachel_data/EndoU_project"
norm_tabs <- fs::path_join(c(basedir, "all_RNA_bg/exp_2"))

all_reads_exp2 <- sum_all_reads(filenames)

usethis::use_data(all_reads_exp2, compress = "xz", overwrite = TRUE)
