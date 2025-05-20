#!/usr/local/bin/Rscript
#
# Jonas Grossmann <jg@fgcz.ethz.ch>
# 2025
#
#


library(stringr)
library(readr)
library(dplyr)
library(tidyr)

desc <- read_delim("DescriptionLinesFromGenomicsFasta.txt", delim = "\t", col_names = FALSE)

desc$X1

descFasta <- data.frame(matrix(ncol = 0, nrow = nrow(desc)))
descFasta$proteinID <- sapply(strsplit(gsub(x = sapply(strsplit(desc$X1, split = " "), function(x)x[1]), pattern = ">lcl\\|CP006706.1_cds_", replacement = ""),split = "_"), function(x)x[1])
descFasta$LocusTag <- gsub(x = desc$X1, pattern = ".*\\[locus_tag=(.*)\\] \\[protein=.*", replacement = "\\1")
descFasta$location <- gsub(x = desc$X1, pattern = ".*\\[location=(.*)\\] \\[gbkey=.*", replacement = "\\1")
descFasta$proteinID2 <- gsub(x = desc$X1, pattern = ".*\\[protein_id=(.*)\\] \\[location=.*", replacement = "\\1")
descFasta$description <- gsub(x = desc$X1, , pattern = ".*\\[protein=(.*)\\] \\[protein_id=.*", replacement = "\\1")

# save as RDS
saveRDS(descFasta, file = "descFasta.rds")
