# 2025-05-19
# this script is used to join the prX data with the methylation data
# storing the results in the resources folder is too big

# load libraries
library(dplyr)
library(readr)
library(tidyr)
library(dplyr)
library(IRanges)


# read in the data
prX <- readRDS("SA6850_prXallWide_moreMetaInfo.rds")
mtX <- readRDS("methylation_data_EGGnogAnnotated.rds")

# globals
IPDrThreshold <- 4

# join mtX 2 prX
prX_wMtX <- left_join(x = prX, y = mtX, by = c("LocusTag" = "locus_tag"))

# write to rds they are too big after the join
#write_rds(x = prX_wMtX, file = "SA6850_prXallWide_moreMetaInfo_methylation.rds")
#write_tsv(prX_wMtX, "SA6850_prXallWide_moreMetaInfo_methylation.tsv")

# join PrX 2 mtx
mtX_wPrX <- left_join(x = mtX, y = prX, by = c("locus_tag" = "LocusTag"))

# GOAL visualize with densities on the chromosome where there are methylations
# maybe use iPDratio as weights

# biocyc table
biocycSA6850 <- read_tsv("Biocyc_allGenesSA6850.txt")


colnames(mtX_wPrX)
unique(mtX_wPrX$group) # use 6850
unique(mtX_wPrX$treatment) # use TSB
# check what is useful in my data
# limit to only one strain (ancestor TSB)
mtX_wPrX_6850_TSB <- mtX_wPrX %>%
  filter(group == "6850") %>%
  filter(treatment == "TSB") %>%
    filter(IPDRatio > IPDrThreshold)


#mtX_wPrX_6850_TSB %>% select(locus_tag, start, strand, feature, mean_iBAQ, AGU, category) |> distinct()
methylation_data <- mtX_wPrX_6850_TSB %>% select(start, strand, feature) |> distinct()

# visualize data as density on the chromosome with start as position
# Load necessary library
library(ggplot2)

# Assuming your data frame is called 'methylation_data'
# Create the density plot with proper direction based on strand
ggplot(methylation_data, aes(x = start)) +
    geom_density(data = subset(methylation_data, strand == "+"),
                 aes(y = after_stat(density), fill = strand),
                 alpha = 0.5) +
    geom_density(data = subset(methylation_data, strand == "-"),
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


# as hist
# Histogram approach
ggplot() +
    geom_histogram(data = subset(methylation_data, strand == "+"),
                   aes(x = start, y = after_stat(count), fill = strand),
                   alpha = 0.5, bins = 500) +
    geom_histogram(data = subset(methylation_data, strand == "-"),
                   aes(x = start, y = -after_stat(count), fill = strand),
                   alpha = 0.5, bins = 500) +
    scale_fill_manual(values = c("+" = "blue", "-" = "red")) +
    labs(
        title = "Distribution of Methylation Sites by Strand (Ancestor TSB)",
        x = "Chromosome Position",
        y = "Count (+ strand up, - strand down)",
        fill = "Strand"
    ) +
    theme_minimal() +
    geom_hline(yintercept = 0, linetype = "dashed")
# save the plot
ggsave("methylation_density_plot.png", width = 10, height = 6)


methylation_data_6850_PASN <- mtX_wPrX %>%
    filter(group == "6850") %>%
    filter(treatment == "PASN") %>%
    filter(IPDRatio > IPDrThreshold) %>%
    select(start, strand, feature) |> distinct()

# Histogram approach
ggplot() +
    geom_histogram(data = subset(methylation_data_6850_PASN, strand == "+"),
                   aes(x = start, y = after_stat(count), fill = strand),
                   alpha = 0.5, bins = 500) +
    geom_histogram(data = subset(methylation_data_6850_PASN, strand == "-"),
                   aes(x = start, y = -after_stat(count), fill = strand),
                   alpha = 0.5, bins = 500) +
    scale_fill_manual(values = c("+" = "blue", "-" = "red")) +
    labs(
        title = "Distribution of Methylation Sites by Strand (Ancestor PASN)",
        x = "Chromosome Position",
        y = "Count (+ strand up, - strand down)",
        fill = "Strand"
    ) +
    theme_minimal() +
    geom_hline(yintercept = 0, linetype = "dashed")
# save the plot
ggsave("methylation_density_plot_6850_PASN.png", width = 10, height = 6)


# how about regions of hypomethylation
head(methylation_data_6850_PASN)

# calc and add column with difference to next position on the same strand
methylation_data_6850_PASN <- methylation_data_6850_PASN %>%
    arrange(start) %>%
    group_by(strand) %>%
    mutate(diff = start - lag(start)) %>%
    ungroup()

# visualize diff per strand as boxplot and use the same colors for strand
ggplot(methylation_data_6850_PASN, aes(x = strand, y = diff, fill = strand)) +
    geom_boxplot() +
    scale_fill_manual(values = c("+" = "blue", "-" = "red")) +
    labs(
        title = "Difference in Chromosome Position by Strand (Ancestor PASN)",
        x = "Strand",
        y = "Difference in Chromosome Position",
        fill = "Strand"
    ) +
    theme_minimal()
# save the plot
ggsave("methylation_diff_boxplot_6850_PASN.png", width = 10, height = 6)


# find regions of hypomethylation (largest diff)
# find the largest difference
largest_diff_plusStrand <- methylation_data_6850_PASN %>%
    filter(diff > 0, strand == "+") %>%
    arrange(desc(diff)) %>%
    mutate(hypoMethRegion_end = start + diff) %>%
    dplyr::slice(1:10) # top 10 largest differences


colnames(biocycSA6850) <- c("geneName", "Acc1", "start", "end", "product", "strand", "GOBP", "GOMF", "pw", "Acc2")

biocycSA6850_slim <- biocycSA6850 %>%
    select(Acc2, start, end, strand) %>%
    mutate(start = as.numeric(start), end = as.numeric(end))


head(biocycSA6850_slim)
head(largest_diff_plusStrand)
# find all genes in the hypomethRegio from biocycSA6850 table

# Function to check if genes overlap with hypo-methylated regions
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

# Filter for just the plus strand genes if your hypometh regions are strand-specific
plus_strand_genes <- biocycSA6850_slim %>%
    filter(strand == "+")

# Find genes overlapping with hypo-methylated regions
overlapping_acc2 <- find_overlapping_genes(plus_strand_genes, largest_diff_plusStrand)

# Print results
cat("Number of Acc2 genes overlapping with hypo-methylated regions:", length(overlapping_acc2), "\n")
print(overlapping_acc2)

# Find the full details of these genes
overlapping_genes_details <- biocycSA6850_slim %>%
    filter(Acc2 %in% overlapping_acc2)

# Print the details
print(overlapping_genes_details)
sortedHypRegionGenes <- overlapping_genes_details[order(overlapping_genes_details$Acc2),]


# this looks more interesting.. here a hypomethylated region.
# https://biocyc.org/tu?orgid=GCF_000462955&id=TU2BZV-39
# RSAU_000070
# RSAU_000071
# RSAU_000072
# RSAU_000073
# RSAU_000074
# RSAU_000075
# RSAU_000076
# RSAU_000077
# RSAU_000078

# RSAU_000253
# RSAU_000254
# RSAU_000255
# RSAU_000256
# RSAU_000257
# RSAU_000258
# RSAU_000259
# RSAU_000262
# RSAU_000264
# RSAU_000266


# for these genes I want to join in the proteomics data
# sortedHypRegionGenes

prXdata <- read_tsv("SA6850_prXallWide_moreMeta.tsv")

POIs <- left_join(x = sortedHypRegionGenes, y = prXdata, by = c("Acc2" = "LocusTag"))

write_tsv(POIs, "hypomethylated_genes_with_prXdata.tsv")
# stringEnrichment
# Now it starts to get interesting:
# https://version-12-0.string-db.org/cgi/network?networkId=bYUQmgbHc0w9




# find regions of hypomethylation (largest diff)
# find the largest difference
largest_diff_AllStrand_all_topX <- methylation_data_6850_PASN %>%
    arrange(desc(diff)) %>%
    mutate(hypoMethRegion_end = start + diff) %>%
    dplyr::slice(1:20) # top 20 largest differences



head(biocycSA6850_slim)
head(largest_diff_AllStrand_all_topX)

# Filter for just the plus strand genes if your hypometh regions are strand-specific
# plus_strand_genes <- biocycSA6850_slim %>%
#     filter(strand == "+")

# Find genes overlapping with hypo-methylated regions
overlapping_acc2_wAll <- find_overlapping_genes(biocycSA6850_slim, largest_diff_AllStrand_all_topX)

# Print results
cat("Number of Acc2 genes overlapping with hypo-methylated regions:", length(overlapping_acc2_wAll), "\n")
print(overlapping_acc2_wAll)

# Find the full details of these genes
overlapping_genes_details_all <- biocycSA6850_slim %>%
    filter(Acc2 %in% overlapping_acc2_wAll)

# Print the details
print(overlapping_genes_details_all)
sortedHypRegionGenes_all <- overlapping_genes_details_all[order(overlapping_genes_details_all$Acc2),]

# for these genes I want to join in the proteomics data
# sortedHypRegionGenes

POIs_all <- left_join(x = sortedHypRegionGenes_all, y = prXdata, by = c("Acc2" = "LocusTag"))

write_tsv(POIs_all, "hypomethylated_genes_with_prXdata_withAll.tsv")
# stringEnrichment
# Now it starts to get interesting:
# redone with all and split for strand
# plus strand
# https://string-db.org/cgi/network?taskId=b8Ce32uVbmpT&sessionId=byObJaqP47cs

# minus strand
# https://string-db.org/cgi/network?taskId=bpsQVJiPeB8s&sessionId=bOEjIAIqtXCK





