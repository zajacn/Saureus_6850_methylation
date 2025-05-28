# 2025-05-28
# this script is used to join the prX data with the methylation data
# storing the results in the resources folder is too big

# load libraries
library(dplyr)
library(readr)
library(tidyr)

# read in the data
prX <- readRDS("SA6850_prXallWide_moreMetaInfo.rds")
mtX <- readRDS("methylation_data_EGGnogAnnotated.rds")


# start here if you are interested in one type of category
# join mtX 2 prX
prX_wMtX <- left_join(x = prX, y = mtX, by = c("LocusTag" = "locus_tag"))



# we know that m6A has usually IPDratios > 3.5
# 5mC has usually IPDratios < 3
# we want to keep all
plot(x = prX_wMtX$IPDRatio, y = log10(prX_wMtX$mean_iBAQ),
     xlab = "IPD ratio", ylab = "mean_iBAQ",
     main = "IPD ratio vs mean_iBAQ",
     col = as.factor(prX_wMtX$feature),
     pch = ".")

colnames(prX_wMtX)

# we do not fully understand idQV yet but according to NZ it is all filtered already
# modified base do not have a value for identificationQv
plot(x = prX_wMtX$IPDRatio, y = prX_wMtX$identificationQv,
     xlab = "IPD ratio", ylab = "identificationQv",
     main = "IPD ratio vs identificationQv",
     col = as.factor(prX_wMtX$feature),
     pch = ".")

# what catetories are there
unique(prX_wMtX$category)

# category of Interest
COI <- "downstream"

table(prX_wMtX$category)
table(prX_wMtX$feature)

# we want to filter out modified_base
prX_wMtX <- prX_wMtX[prX_wMtX$feature != "modified_base",]
table(prX_wMtX$feature)

# filter for Ancestor only
unique(prX_wMtX$group)

# anyway we are only looking into ancestor
prX_wMtX <- prX_wMtX[prX_wMtX$group == "6850",]

# now summarize per AGU and treatment
table(prX_wMtX$feature)
table(prX_wMtX$category)

# we are only interested in upstream features
prX_wMtX <- prX_wMtX[prX_wMtX$category == COI,]
table(prX_wMtX$feature)
table(prX_wMtX$category)

# stringNW and enrichment:
# annotated proteome: STRG0A52IKN

# prX_wMtX %>% reframe(AGU, max_IPDratio = max(IPDRatio, na.rm = TRUE), nCount = n())

# Summarize data to get max IPDRatio per AGU and treatment
countAll_perAGUnTreatment <- prX_wMtX %>%
    group_by(AGU, treatment) %>%
    summarize(nCount = n(), .groups = "drop")


# Reshape to wide format (optional)
countAll_perAGUnTreatment_wide <-  prX_wMtX %>%
    group_by(AGU, treatment) %>%
    summarize(nCount = n(), .groups = "drop") %>%
    pivot_wider(names_from = treatment, values_from = nCount)

# some do not have any counts (NA) -> impute small value to not loose in for the log2
countAll_perAGUnTreatment_wide <- countAll_perAGUnTreatment_wide %>%
    mutate(PASN = ifelse(is.na(PASN), 0, PASN),
           TSB = ifelse(is.na(TSB), 0, TSB)) |> select(AGU, PASN, TSB)

# Replace NA with a small value to avoid log2(0) issues
countAll_perAGUnTreatment_wide <- na.omit(countAll_perAGUnTreatment_wide)

# impute -> impute small value to not loose in for the log2 for zeros
countAll_perAGUnTreatment_wide <- countAll_perAGUnTreatment_wide %>%
    mutate(PASN = ifelse(PASN == 0, 0.01, PASN),
           TSB = ifelse(TSB == 0, 0.01, TSB))

# Calculate log2 ratio
count_log2Ratio <- countAll_perAGUnTreatment_wide %>%
    mutate(log2_ratio = log2(PASN / TSB)) |> select(AGU, log2_ratio)

hist(count_log2Ratio$log2_ratio,
     xlab = "log2 Ratio of PASN/TSB", ylab = "Frequency",
     main = "Distribution of log2 Ratio of PASN/TSB",
     col = "lightblue", border = "black", breaks=50)

# try to look into GSEA with all types and features in
#COI <- "all"
(fn <- paste("GSEA_log2CountPASNvsTSB_",COI,"_features", ".txt", sep = ""))

write.table(count_log2Ratio, file = fn, sep = "\t", row.names = FALSE, quote = FALSE)

# join in proteomics data
colnames(prX)
PrXforLog2 <- prX %>%
     select(AGU, protein_length, nrPeptides, mean_iBAQ, statistic.PASNvsTSB_givenAncestor)

# some plots
# some corr
plot(x = PrXforLog2$protein_length, y = PrXforLog2$nrPeptides,
     xlab = "protLenght", ylab = "nrPeptides",
     main = "Length vs nrPeptpides",
     pch = ".")

# corr not there as expected
plot(x = sqrt(PrXforLog2$protein_length), y = log10(PrXforLog2$mean_iBAQ),
     xlab = "protLenght", ylab = "iBAQ",
     main = "Length vs nrPeptpides",
     pch = ".")

# join in log2 ratio
PrXforLog2 <- left_join(x = PrXforLog2, y = count_log2Ratio, by = "AGU")

# Do we see structure
plot(x = PrXforLog2$log2_ratio, y = PrXforLog2$statistic.PASNvsTSB_givenAncestor,
     xlab = "log2 Ratio of PASN/TSB", ylab = "statistic PASN vs TSB",
     main = "log2 Ratio vs statistic PASN vs TSB",
     pch = "x")
abline(h=0, col = "red", lty = 2)
abline(v=0, col = "green", lty = 2)


# disappointingly no structure for all the feature types ;(
# how about the proteins that do NOT have a log2Ratio but we see them expressed
sum(is.na(PrXforLog2$log2_ratio))
plot(x = PrXforLog2$mean_iBAQ[is.na(PrXforLog2$log2_ratio)], y = PrXforLog2$statistic.PASNvsTSB_givenAncestor[is.na(PrXforLog2$log2_ratio)],
     xlab = "log2 Ratio of PASN/TSB", ylab = "statistic PASN vs TSB",
     main = "log2 Ratio vs statistic PASN vs TSB",
     pch = "x")

# AGU expressed but not methylated
write.table(PrXforLog2$AGU[is.na(PrXforLog2$log2_ratio)],
            file = "AGU_expressed_but_not_methylated.txt",
            sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)




