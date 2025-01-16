#!/usr/local/bin/Rscript
#
# Jonas Grossmann <jg@fgcz.ethz.ch>
# 2024
#
#

# Jonas looks into methylation data
library(tidyverse)


mtx <- readRDS("resources/methylation_data.rds")
dim(mtx)


# try to visualize mA6 as binary heatmap along the genome
str(mtx)

mtx |> select(c(treatment, group, strand))


  as("data.frame") |>
  as.matrix() |>
  pheatmap::pheatmap()

