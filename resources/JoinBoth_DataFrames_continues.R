# 2025-04-11
# this script is used to join the prX data with the methylation data
# storing the results in the resources folder is too big

# load libraries
library(dplyr)
library(readr)
library(tidyr)

# read in the data
prX <- readRDS("SA6850_prXallWide_moreMetaInfo.rds")
mtX <- readRDS("methylation_data_EGGnogAnnotated.rds")

# join mtX 2 prX
prX_wMtX <- left_join(x = prX, y = mtX, by = c("LocusTag" = "locus_tag"))

# write to rds they are too big after the join
#write_rds(x = prX_wMtX, file = "SA6850_prXallWide_moreMetaInfo_methylation.rds")
#write_tsv(prX_wMtX, "SA6850_prXallWide_moreMetaInfo_methylation.tsv")

# join PrX 2 mtx
mtX_wPrX <- left_join(x = mtX, y = prX, by = c("locus_tag" = "LocusTag"))



# look into aspect of avg protein expression and methylation

prX_wMtX$IPDRatio
prX_wMtX$mean_iBAQ

# IPD ratio and iBAQ
colnames(prX_wMtX)
plot(x = prX_wMtX$IPDRatio, y = log10(prX_wMtX$mean_iBAQ),
     xlab = "IPD ratio", ylab = "mean_iBAQ",
     main = "IPD ratio vs mean_iBAQ",
     pch = ".")

# smoothscatter
smoothScatter(x = prX_wMtX$IPDRatio, y = log10(prX_wMtX$mean_iBAQ),
               xlab = "IPD ratio", ylab = "mean_iBAQ",
               main = "IPD ratio vs mean_iBAQ",
               pch = ".")


# Interesting pattern visible that some expressed proteins have IPDratios >4
# What are these? any functional enrichment?

IPDrThreshold <- 4
Methylated <- prX_wMtX[prX_wMtX$IPDRatio > IPDrThreshold,]

methylatedProteins <- unique(Methylated$proteinID)
write.table(methylatedProteins, file = "methylatedProteins.txt",
            sep = "\t", row.names = FALSE, col.names = FALSE,
            quote = FALSE)

# further split methylated proteins for
table(Methylated$feature)

m4C_proteins <- unique(Methylated$proteinID[Methylated$feature == "m4C"])
m6A_proteins <- unique(Methylated$proteinID[Methylated$feature == "m6A"])

table(Methylated$category)

length(na.omit(unique(Methylated$proteinID[Methylated$category == "upstream"])))
write.table(na.omit(unique(Methylated$proteinID[Methylated$category == "upstream"])),
            file = "methylatedProteins_upstream.txt",
            sep = "\t", row.names = FALSE, col.names = FALSE,
            quote = FALSE)
# https://version-12-0.string-db.org/cgi/network?networkId=bnfYns5X67qY

length(na.omit(unique(Methylated$proteinID[Methylated$category == "gene"])))
write.table(na.omit(unique(Methylated$proteinID[Methylated$category == "gene"])),
            file = "methylatedProteins_gene.txt",
            sep = "\t", row.names = FALSE, col.names = FALSE,
            quote = FALSE)
# https://version-12-0.string-db.org/cgi/network?networkId=bF5adZzioD3B


length(na.omit(unique(Methylated$proteinID[Methylated$category == "downstream"])))
write.table(na.omit(unique(Methylated$proteinID[Methylated$category == "downstream"])),
            file = "methylatedProteins_downstream.txt",
            sep = "\t", row.names = FALSE, col.names = FALSE,
            quote = FALSE)
# https://version-12-0.string-db.org/cgi/network?networkId=bmIktug2JDNu


# stringNW and enrichment:
# annotated proteome: STRG0A52IKN
# https://version-12-0.string-db.org/cgi/network?networkId=bhZ55QMAsfUC
# w/ adjusted BG
# https://version-12-0.string-db.org/cgi/network?networkId=bbs5Y269mlOr


# Look into "most methylated proteins".. GSEA ranked for IPD ratio and split for upstream, gene, downstream

# upstream methylated
upstream_methylated <- prX_wMtX[prX_wMtX$IPDRatio > IPDrThreshold &
                                    prX_wMtX$category == "upstream",]
# max_IPDratio <- prX_wMtX %>%
#group_by(proteinID) %>%
#    summarise(max_IPDratio = max(IPDRatio, na.rm = TRUE))

upstream_methylated %>% reframe(AGU, max_IPDratio = max(IPDRatio, na.rm = TRUE), nCount = n()) |>
    distinct()

# Summarize data to get max IPDRatio per AGU and treatment
max_ipd_perAGUnTreatment <- upstream_methylated %>%
    group_by(AGU, treatment) %>%
    summarize(max_IPDRatio = max(IPDRatio), .groups = "drop")

# I would like to generate a log2 Ratio for maxIPDratio of the different treatments (PASN_vs_TSB)
# Reshape to wide format (optional)
max_ipd_wide <- max_ipd_perAGUnTreatment %>%
    pivot_wider(names_from = treatment, values_from = max_IPDRatio)

# Calculate log2 ratio
max_ipd_log2Ratio <- max_ipd_wide %>%
    mutate(log2_ratio = log2(PASN / TSB)) |> select(AGU, log2_ratio)
write.table(max_ipd_log2Ratio, file = "max_ipd_log2Ratio.txt",
            sep = "\t", row.names = FALSE, col.names = TRUE,
            quote = FALSE)



# within gene methylated
withinGene_methylated <- prX_wMtX[prX_wMtX$IPDRatio > IPDrThreshold &
                                    prX_wMtX$category == "gene",]

# Summarize data to get max IPDRatio per AGU and treatment
max_ipd_perAGUnTreatment_gene <- withinGene_methylated %>%
    group_by(AGU, treatment) %>%
    summarize(max_IPDRatio = max(IPDRatio), .groups = "drop")

# I would like to generate a log2 Ratio for maxIPDratio of the different treatments (PASN_vs_TSB)
# Reshape to wide format (optional)
max_ipd_wide_gene <- max_ipd_perAGUnTreatment_gene %>%
    pivot_wider(names_from = treatment, values_from = max_IPDRatio)

# Calculate log2 ratio
max_ipd_log2Ratio_gene <- max_ipd_wide_gene %>%
    mutate(log2_ratio = log2(PASN / TSB)) |> select(AGU, log2_ratio)
write.table(max_ipd_log2Ratio_gene, file = "max_ipd_log2Ratio_gene.txt",
            sep = "\t", row.names = FALSE, col.names = TRUE,
            quote = FALSE)

# write out some more

AGUmethInPASN <- max_ipd_wide_gene |> select(AGU, PASN)
AGUmethInTSB <- max_ipd_wide_gene |> select(AGU, TSB)

write.table(AGUmethInPASN, file = "AGUmethInPASN.txt",
            sep = "\t", row.names = FALSE, col.names = TRUE,
            quote = FALSE)
write.table(AGUmethInTSB, file = "AGUmethInTSB.txt",
            sep = "\t", row.names = FALSE, col.names = TRUE,
            quote = FALSE)

# no string enrichment
# https://version-12-0.string-db.org/cgi/globalenrichment?networkId=bC30fyhuXxZB

# what are the AGU with largest maxIPDratio for PASN
max_ipd_wide_gene %>% select(AGU, PASN) %>%
    arrange(desc(PASN)) %>%
    slice_head(n = 50) %>%
    write.table(file = "topX_max_ipd_PASN.txt",
                sep = "\t", row.names = FALSE, col.names = TRUE,
                quote = FALSE)

# what are the AGU with largest maxIPDratio for PASN
max_ipd_wide_gene %>% select(AGU, TSB) %>%
    arrange(desc(TSB)) %>%
    slice_head(n = 50) %>%
    write.table(file = "topX_max_ipd_TSB.txt",
                sep = "\t", row.names = FALSE, col.names = TRUE,
                quote = FALSE)

# no enrichements but:
# https://version-12-0.string-db.org/cgi/network?networkId=bwxieg8S7WEs
# unknown network (max IPD for PASN)

# no enrichment but..
# https://version-12-0.string-db.org/cgi/network?networkId=blSaz901aEBe

# summarize how often a protein is methylated
Topcount_Methylated_gene <- prX_wMtX[prX_wMtX$IPDRatio > IPDrThreshold &
             prX_wMtX$category == "gene",] %>%
    select(proteinID, AGU, IPDRatio) %>% group_by(AGU) %>%
    summarize(count = n(), .groups = "drop") %>% arrange(desc(count))
write.table(Topcount_Methylated_gene, file = "Topcount_Methylated.txt",
            sep = "\t", row.names = FALSE, col.names = TRUE,
            quote = FALSE)

table(prX_wMtX$category)
# summarize how often a protein is methylated
Topcount_Methylated_upstream <- prX_wMtX[prX_wMtX$IPDRatio > IPDrThreshold &
                                    prX_wMtX$category == "upstream",] %>%
    select(proteinID, AGU, IPDRatio) %>% group_by(AGU) %>%
    summarize(count = n(), .groups = "drop") %>% arrange(desc(count))
write.table(Topcount_Methylated_upstream, file = "Topcount_Methylated_upstream.txt",
            sep = "\t", row.names = FALSE, col.names = TRUE,
            quote = FALSE)

# summarize how often a protein is methylated
Topcount_Methylated_dwnstream <- prX_wMtX[prX_wMtX$IPDRatio > IPDrThreshold &
                                             prX_wMtX$category == "downstream",] %>%
    select(proteinID, AGU, IPDRatio) %>% group_by(AGU) %>%
    summarize(count = n(), .groups = "drop") %>% arrange(desc(count))
write.table(Topcount_Methylated_dwnstream, file = "Topcount_Methylated_dwnstream.txt",
            sep = "\t", row.names = FALSE, col.names = TRUE,
            quote = FALSE)



# how about the other way around?
nonmethylatedProteins <- unique(prX_wMtX$proteinID[prX_wMtX$IPDRatio < IPDrThreshold])
write.table(nonmethylatedProteins, file = "nonmethylatedProteins.txt",
            sep = "\t", row.names = FALSE, col.names = FALSE,
            quote = FALSE)

# ORAs
# let's split the expressed and methylated proteins in several groups according to iBAQ

# background
allIdentifiedProteins <- unique(prX_wMtX$proteinID)
write.table(allIdentifiedProteins, file = "allIdentifiedProteins.txt",
            sep = "\t", row.names = FALSE, col.names = FALSE,
            quote = FALSE)


# find useful thresholds based on expression
# get the distribution of iBAQ values and split in 5 equal groups

# get the distribution of iBAQ values
hist(log10(prX_wMtX$mean_iBAQ), breaks = 100,
     xlab = "log10(mean iBAQ)", ylab = "Frequency",
     main = "Distribution of iBAQ values",
     col = "lightblue", border = "black")

# get the quantiles for log10
quantiles_log10 <- quantile(log10(prX_wMtX$mean_iBAQ), probs = seq(0, 1, 0.2))
quantiles_log10

# split the data into 5 groups
prX_wMtX$group <- cut(log10(prX_wMtX$mean_iBAQ),
                       breaks = quantiles_log10,
                       labels = c("very low", "low", "medium", "high", "very high"),
                       include.lowest = TRUE)

table(prX_wMtX$group)
# get the methylated proteins in each group
methylated_proteins_group <- prX_wMtX[prX_wMtX$IPDRatio > IPDrThreshold,]
methylated_proteins_group$group <- cut(log10(methylated_proteins_group$mean_iBAQ),
                                         breaks = quantiles_log10,
                                         labels = c("very low", "low", "medium", "high", "very high"),
                                         include.lowest = TRUE)
table(methylated_proteins_group$group)

# write the methylated proteins in each group to a file
length(unique(methylated_proteins_group$proteinID[methylated_proteins_group$group == "very low"]))
write.table(unique(methylated_proteins_group$proteinID[methylated_proteins_group$group == "very low"]),
            file = "methylated_proteins_very_low.txt",
            sep = "\t", row.names = FALSE, col.names = FALSE,
            quote = FALSE)
# https://version-12-0.string-db.org/cgi/network?networkId=bLyiGOUWYtZ3

length(unique(methylated_proteins_group$proteinID[methylated_proteins_group$group == "low"]))
write.table(unique(methylated_proteins_group$proteinID[methylated_proteins_group$group == "low"]),
            file = "methylated_proteins_low.txt",
            sep = "\t", row.names = FALSE, col.names = FALSE,
            quote = FALSE)

length(unique(methylated_proteins_group$proteinID[methylated_proteins_group$group == "medium"]))
write.table(unique(methylated_proteins_group$proteinID[methylated_proteins_group$group == "medium"]),
            file = "methylated_proteins_medium.txt",
            sep = "\t", row.names = FALSE, col.names = FALSE,
            quote = FALSE)
# https://version-12-0.string-db.org/cgi/network?networkId=b1IakohUdFBC


length(unique(methylated_proteins_group$proteinID[methylated_proteins_group$group == "high"]))
write.table(unique(methylated_proteins_group$proteinID[methylated_proteins_group$group == "high"]),
            file = "methylated_proteins_high.txt",
            sep = "\t", row.names = FALSE, col.names = FALSE,
            quote = FALSE)

length(unique(methylated_proteins_group$proteinID[methylated_proteins_group$group == "very high"]))
write.table(unique(methylated_proteins_group$proteinID[methylated_proteins_group$group == "very high"]),
            file = "methylated_proteins_very_high.txt",
            sep = "\t", row.names = FALSE, col.names = FALSE,
            quote = FALSE)
# https://version-12-0.string-db.org/cgi/network?networkId=bcIvwiSQp4sc

# for doing GSEA get mean_IPDratio for each proteinID
# get the mean IPD ratio for each proteinID
mean_IPDratio <- prX_wMtX %>%
  group_by(proteinID) %>%
  summarise(mean_IPDratio = mean(IPDRatio, na.rm = TRUE))
write.table(mean_IPDratio, file = "mean_IPDratio.txt",
            sep = "\t", row.names = FALSE, col.names = TRUE,
            quote = FALSE)

# get the max IPD ratio for each proteinID
max_IPDratio <- prX_wMtX %>%
  group_by(proteinID) %>%
  summarise(max_IPDratio = max(IPDRatio, na.rm = TRUE))
write.table(na.omit(max_IPDratio), file = "max_IPDratio.txt",
            sep = "\t", row.names = FALSE, col.names = TRUE,
            quote = FALSE)

# summarize IPD ratio for each proteinID how often is it >4
IPDratio_summary <- prX_wMtX %>% select(proteinID, IPDRatio) %>% distinct() %>%
  group_by(proteinID) %>%
  summarise(IPDratio_count = sum(IPDRatio > IPDrThreshold, na.rm = TRUE))
write.table(na.omit(IPDratio_summary), file = "IPDratio_summary.txt",
            sep = "\t", row.names = FALSE, col.names = TRUE,
            quote = FALSE)

# nothing in GSEA
# https://version-12-0.string-db.org/cgi/globalenrichment?networkId=bvIRPq9Bxq4h


# further exploration
plot(prX_wMtX$identificationQv ~ prX_wMtX$IPDRatio)
plot(prX_wMtX$identificationQv ~ log10(prX_wMtX$mean_iBAQ))

plot(prX_wMtX$identificationQv ~ prX_wMtX$IPDRatio,
     xlab = "IPD ratio", ylab = "identification Qv",
     main = "IPD ratio vs identification Qv",
     pch = ".", cex = (log10(prX_wMtX$mean_iBAQ+0.1)/2), col = "blue")




# try AI
# https://app.formulabot.com/

# reduce df to only ancestor
unique(prX_wMtX$Name)

prX_wMtX_ancestor <- prX_wMtX[grepl(x = prX_wMtX$Name, pattern = "6850_"),]
write_tsv(x = prX_wMtX_ancestor, file = "SA6850_prXallWide_moreMetaInfo_methylation_ancestor.tsv")



