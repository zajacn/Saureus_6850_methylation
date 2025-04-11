# 2025-04-11
# this script is used to join the prX data with the methylation data
# storing the results in the resources folder is too big

# read in the data
prX <- readRDS("resources/SA6850_prXallWide_moreMetaInfo.rds")
mtX <- readRDS("resources/methylation_data_EGGnogAnnotated.rds")

# join mtX 2 prX
library(dplyr)
prX_wMtX <- left_join(x = prX, y = mtX, by = c("LocusTag" = "locus_tag"))

# write to rds
write_rds(x = prX_wMtX, file = "resources/SA6850_prXallWide_moreMetaInfo_methylation.rds")
write_tsv(prX_wMtX, "resources/SA6850_prXallWide_moreMetaInfo_methylation.tsv")

# join PrX 2 mtx
mtX_wPrX <- left_join(x = mtX, y = prX, by = c("locus_tag" = "LocusTag"))



