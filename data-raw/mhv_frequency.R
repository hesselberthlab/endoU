library(endoU)
library(tidyverse)
library(fs)
library(usethis)


basedir <- "~/Dropbox (Hesselberth Lab)/Rachel_data/EndoU_project"
dinuc <- c("AA","GA", "CA", "UA","AC", "CC", "GC", "UC", "AG", "CG", "GG", "UG", "AU", "CU", "GU", "UU")
df_dinuc <- data_frame(dinuc)

background_dinuc_counts <- fs::path_join(c(basedir, "mhv_fasta_dinuc_counts.txt"))

mhv_freq <- readr::read_tsv(background_dinuc_counts) %>%
  dplyr::mutate(frequency = count/sum(count)) %>%
  mutate(count_sum = sum(count))

usethis::use_data(mhv_freq, compress = "xz", overwrite = TRUE)


