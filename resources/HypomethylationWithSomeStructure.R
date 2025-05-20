# 2025-05-19
# this script is used to join the prX data with the methylation data
# storing the results in the resources folder is too big

# load libraries
library(dplyr)
library(readr)
library(tidyr)
library(dplyr)
library(IRanges)
library(ggplot2)


# Function to check if genes overlap with hypo-methylated regions
# this function does not care about the direction or coding strand
find_overlapping_genes <- function(genes_df, hypometh_df) {
    # Create ranges for genes
    gene_ranges <- IRanges(start = genes_df$start,
                           end = genes_df$end,
                           names = genes_df$Acc2)

    # Create ranges for hypo-methylated regions
    # Assuming start column is the beginning of the region
    hypometh_ranges <- IRanges(start = hypometh_df$start,
                               end = hypometh_df$hypoMethRegion_end)

    # Find overlaps
    overlaps <- findOverlaps(gene_ranges, hypometh_ranges)

    # Get the names of overlapping genes
    overlapping_genes <- names(gene_ranges)[queryHits(overlaps)]

    # Return unique gene names
    return(unique(overlapping_genes))
}


#
# data to read in
#


# read in the data
prX <- readRDS("SA6850_prXallWide_moreMetaInfo.rds")
mtX <- readRDS("methylation_data_EGGnogAnnotated.rds")
# biocyc table here all genes are listed from left to rigth - no matter if + or - strand
biocycSA6850 <- read_tsv("Biocyc_allGenesSA6850.txt")
colnames(biocycSA6850) <- c("geneName", "Acc1", "start", "end", "product", "strand", "GOBP", "GOMF", "pw", "Acc2")

biocycSA6850_slim <- biocycSA6850 %>%
    select(Acc2, start, end, strand) %>%
    mutate(start = as.numeric(start), end = as.numeric(end))



# globals
IPDrThreshold <- 4
SpacingCutoff <- 25000

# join mtX 2 prX
# prX_wMtX <- left_join(x = prX, y = mtX, by = c("LocusTag" = "locus_tag"))
mtX_wPrX <- left_join(x = mtX, y = prX, by = c("locus_tag" = "LocusTag"))

# GOAL visualize with densities on the chromosome where there are methylations
colnames(mtX_wPrX)
unique(mtX_wPrX$group) # use 6850
unique(mtX_wPrX$treatment) # use TSB

# check what is useful in my data
# limit to only one strain (ancestor TSB)
clone_of_interest <- "6850"
condition_of_interest <- "PASN"

# filter the data for the specific clone and condition
mtX_wPrX_OfInterest <- mtX_wPrX %>%
    filter(group == clone_of_interest) %>%
    filter(treatment == condition_of_interest) %>%
    filter(IPDRatio > IPDrThreshold)


# slim
methylation_slim <- mtX_wPrX_OfInterest %>% select(start, strand, feature) |> distinct()

# visualize data as density on the chromosome with start as position
# Create the density plot with proper direction based on strand
ggplot(methylation_slim, aes(x = start)) +
    geom_density(data = subset(methylation_slim, strand == "+"),
                 aes(y = after_stat(density), fill = strand),
                 alpha = 0.5) +
    geom_density(data = subset(methylation_slim, strand == "-"),
                 aes(y = -after_stat(density), fill = strand),
                 alpha = 0.5) +
    scale_fill_manual(values = c("+" = "blue", "-" = "red")) +
    labs(
        title = "Density of m6A Methylation Sites by Strand",
        x = "Chromosome Position",
        y = "Density",
        fill = "Strand"
    ) +
    theme_minimal() +
    geom_hline(yintercept = 0, linetype = "dashed")


# as hist -> much better
# Histogram approach
ggplot() +
    geom_histogram(data = subset(methylation_slim, strand == "+"),
                   aes(x = start, y = after_stat(count), fill = strand),
                   alpha = 0.5, bins = 500) +
    geom_histogram(data = subset(methylation_slim, strand == "-"),
                   aes(x = start, y = -after_stat(count), fill = strand),
                   alpha = 0.5, bins = 500) +
    scale_fill_manual(values = c("+" = "blue", "-" = "red")) +
    labs(
        title = paste0("Distribution of Methylation Sites by Strand (",clone_of_interest," ",condition_of_interest,")"),
        x = "Chromosome Position",
        y = "Count (+ strand up, - strand down)",
        fill = "Strand"
    ) +
    theme_minimal() +
    geom_hline(yintercept = 0, linetype = "dashed")

(fN <- paste("Methylation_hist_on_Chromosome_", clone_of_interest, "_", condition_of_interest, ".png", sep = ""))
# save the plot
ggsave(fN, width = 10, height = 6)



#
#
#  Go for regions that are NOT methylated
#
#


# calc and add column with difference to next position on the same strand
# double checked
methylation_data_with_diff <- methylation_slim %>%
    arrange(start) %>%
    group_by(strand) %>%
    mutate(diff = start - lag(start)) %>%
    ungroup()

# visualize diff per strand as boxplot and use the same colors for strand
ggplot(methylation_data_with_diff, aes(x = strand, y = diff, fill = strand)) +
    geom_boxplot() +
    scale_fill_manual(values = c("+" = "blue", "-" = "red")) +
    labs(
        title = paste0("Spacing between methylation sites (",clone_of_interest," ",condition_of_interest,")"),
        x = "Strand",
        y = "Difference in Chromosome Position",
        fill = "Strand"
    ) +
    theme_minimal()

# save the plot
(fN <- paste("Methylation_spacing_on_Chromosome_", clone_of_interest, "_", condition_of_interest, ".png", sep = ""))
ggsave(fN, width = 10, height = 6)


# find regions of hypomethylation (largest diff)
# find the largest difference
methylation_withDiffs <- methylation_data_with_diff %>%
    arrange(desc(diff)) %>%
    mutate(hypoMethRegion_end = start + diff)


# spacing regions of interest
largest_diff_AllStrand_ofInterest <- methylation_withDiffs %>%
    filter(diff > SpacingCutoff) %>%
    select(start, hypoMethRegion_end, diff, strand, feature)

(fN <- paste("Large_Spaces_", clone_of_interest, "_", condition_of_interest, ".tsv", sep = ""))
write_tsv(largest_diff_AllStrand_ofInterest, fN)

# Find genes overlapping with hypo-methylated regions
genes_in_hypo_regions <- find_overlapping_genes(biocycSA6850_slim, largest_diff_AllStrand_ofInterest)

# Some genes in biocyc do not have an RSAU_Number -> but NA --> follow up on these more?? what are these?
genes_in_hypo_regions <- genes_in_hypo_regions[!is.na(genes_in_hypo_regions)]


# Print results
cat("Number of Acc2 genes overlapping with hypo-methylated regions:", length(genes_in_hypo_regions), "\n")
print(genes_in_hypo_regions)

# Find the full details of these genes
overlapping_genes_details_all <- biocycSA6850_slim %>%
    filter(Acc2 %in% genes_in_hypo_regions)

# Print the details
print(overlapping_genes_details_all)

# bring in AGU
desc <- readRDS("descFasta.rds")
overlapping_genes_details_all <- left_join(overlapping_genes_details_all, desc, by = c("Acc2" = "LocusTag"))

# save the results
(fN <- paste("Hypomethylated_genes_n_proteins_", clone_of_interest, "_", condition_of_interest, ".tsv", sep = ""))
write_tsv(overlapping_genes_details_all, fN)

# also join in proteomics expression data
prXdata <- read_tsv("SA6850_prXallWide_moreMeta.tsv")
genes_w_proteinExpression_in_hypoRegions <- left_join(x = overlapping_genes_details_all, y = prXdata, by = c("Acc2" = "LocusTag"))

(fN <- paste("Hypomethylated_genes_n_proteins_", clone_of_interest, "_", condition_of_interest, "_with_proteomics.tsv", sep = ""))
write_tsv(genes_w_proteinExpression_in_hypoRegions, fN)

# SA6850 -> TSB
# minus strand
# https://version-12-0.string-db.org/cgi/network?networkId=bkbClZMyozX6

# plus strand
# https://version-12-0.string-db.org/cgi/network?networkId=bZclsyBXMdc8


# SA6850 -> PASN
# minus strand
# https://version-12-0.string-db.org/cgi/network?networkId=byFYmnWobHl6

# plus strand
# https://version-12-0.string-db.org/cgi/network?networkId=bZclsyBXMdc8



#
# SB0804 -> PASN
# minus strand
# https://version-12-0.string-db.org/cgi/network?networkId=bYdYaaUmS92A
# plus strand
# https://version-12-0.string-db.org/cgi/network?networkId=bSeB33oVQsQT


# compare my iBAQs in general to iBAQ in hyporegions

#
# Create the density plot with proper direction based on strand
ggplot(genes_w_proteinExpression_in_hypoRegions, aes(x = start)) +
    geom_density(data = subset(genes_w_proteinExpression_in_hypoRegions, strand == "+"),
                 aes(y = after_stat(density), fill = strand),
                 alpha = 0.5) +
    geom_density(data = subset(genes_w_proteinExpression_in_hypoRegions, strand == "-"),
                 aes(y = -after_stat(density), fill = strand),
                 alpha = 0.5) +
    scale_fill_manual(values = c("+" = "blue", "-" = "red")) +
    labs(
        title = "Density of hypo regions by Strand",
        x = "Chromosome Position",
        y = "Density",
        fill = "Strand"
    ) +
    theme_minimal() +
    geom_hline(yintercept = 0, linetype = "dashed")
fN <- paste("HypoRegions_density_on_Chromosome_", clone_of_interest, "_", condition_of_interest, ".png", sep = "")
# save the plot
ggsave(fN, width = 10, height = 6)


# as hist -> much better
# Histogram approach
ggplot() +
    geom_histogram(data = subset(genes_w_proteinExpression_in_hypoRegions, strand == "+"),
                   aes(x = start, y = after_stat(count), fill = strand),
                   alpha = 0.5, bins = 500) +
    geom_histogram(data = subset(genes_w_proteinExpression_in_hypoRegions, strand == "-"),
                   aes(x = start, y = -after_stat(count), fill = strand),
                   alpha = 0.5, bins = 500) +
    scale_fill_manual(values = c("+" = "blue", "-" = "red")) +
    labs(
        title = paste0("Distribution of Methylation Sites by Strand (",clone_of_interest," ",condition_of_interest,")"),
        x = "Chromosome Position",
        y = "Count (+ strand up, - strand down)",
        fill = "Strand"
    ) +
    theme_minimal() +
    geom_hline(yintercept = 0, linetype = "dashed")
(fN <- paste("HypoRegions_hist_on_Chromosome_", clone_of_interest, "_", condition_of_interest, ".png", sep = ""))
# save the plot
ggsave(fN, width = 10, height = 6)


hist(log10(prXdata$mean_iBAQ))
hist(log10(genes_w_proteinExpression_in_hypoRegions$mean_iBAQ))

# combine these 2 histograms
ggplot() +
    geom_histogram(data = prXdata,
                   aes(x = log10(mean_iBAQ), y = after_stat(count), fill = "All Proteins"),
                   alpha = 0.5, bins = 20) +
    geom_histogram(data = genes_w_proteinExpression_in_hypoRegions,
                   aes(x = log10(mean_iBAQ), y = after_stat(count), fill = "Hypo Methylated Proteins"),
                   alpha = 0.5, bins = 20) +
    scale_fill_manual(values = c("All Proteins" = "blue", "Hypo Methylated Proteins" = "red")) +
    labs(
        title = paste0("Distribution of iBAQ values (log10) (",clone_of_interest," ",condition_of_interest,")"),
        x = "log10(iBAQ)",
        y = "Count",
        fill = "Proteins"
    ) +
    theme_minimal() +
    geom_hline(yintercept = 0, linetype = "dashed")
(fN <- paste("HypoRegions_iBAQ_vsAll_", clone_of_interest, "_", condition_of_interest, ".png", sep = ""))
# save the plot
ggsave(fN, width = 10, height = 6)
