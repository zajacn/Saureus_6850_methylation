#!/usr/local/bin/Rscript
#
# Jonas Grossmann <jg@fgcz.ethz.ch>
# 2025
#
#
# 2025-05-28
# this script is used to join the prX data with the methylation data

# load libraries
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(purrr)

# set working directory to source file location
# read in the data
prX <- readRDS("SA6850_prXallWide_moreMetaInfo.rds")
mtX <- readRDS("methylation_data_EGGnogAnnotated.rds")

# here mtX centric
mtX_wPrX <- left_join(x = mtX, y = prX, by = c("locus_tag" = "LocusTag"))

# what catetories are there
table(mtX_wPrX$category)

# category of Interest
# no filtering here since we wanna look at all methyation marks
COI <- "allCategories" # this will also be the description of the files outputed


# we want to filter out modified_base -> it is unclear what type of feature it is
mtX_wPrX <- mtX_wPrX[mtX_wPrX$feature != "modified_base",]
table(mtX_wPrX$feature)
table(mtX_wPrX$category)

# filter for Ancestor only
unique(mtX_wPrX$group)

# anyway we are only looking into ancestor (the rest is for students of Bio253)
mtX_wPrX <- mtX_wPrX[mtX_wPrX$group == "6850",]

# HERE we are only interested in one type features
# if multiple once are to be kept, filter out unwanted once
#mtX_wPrX <- mtX_wPrX[mtX_wPrX$category == COI,]
table(mtX_wPrX$feature) # what features 6mA and 4mC should be only
table(mtX_wPrX$category) # here what marks are kept to be analyzed and viz

# we are only interested in the TSB treatment
#condition_of_interest <- "TSB"
#mtX_wPrX_OneTreatment <- mtX_wPrX[mtX_wPrX$treatment == condition_of_interest,]


# make df more manageable
# slim
methylation_slim <- mtX_wPrX %>% select(treatment, start, strand, feature, locus_tag) |> distinct() # here we do have multiple reps in

# for the plotting
clone_of_interest <- "Ancestor6850"


smoothin_param <- 0.001
# visualize data as density on the chromosome with start as position
# Create the density plot with proper direction based on strand
ggplot(methylation_slim, aes(x = start)) +
    geom_density(data = subset(methylation_slim, strand == "+"),
                 aes(y = after_stat(density), fill = strand),
                 alpha = 0.5, adjust = smoothin_param) +  # <-- less smoothing) +
    geom_density(data = subset(methylation_slim, strand == "-"),
                 aes(y = -after_stat(density), fill = strand),
                 alpha = 0.5,
                 adjust = smoothin_param) +  # <-- less smoothing) +
    scale_fill_manual(values = c("+" = "blue", "-" = "red")) +
    facet_grid(treatment ~ .) +
    labs(
        title = "Density of Methylation Sites by Strand",
        x = "Chromosome Position",
        y = "Density",
        fill = "Strand"
    ) +
    theme_minimal() +
    geom_hline(yintercept = 0, linetype = "dashed")
(fN <- paste("Methylation_density_on_Chromosome_v2", clone_of_interest, "_category_", COI, ".png", sep = ""))
# save the plot
ggsave(fN, width = 20, height = 6)

# as hist -> much better
# Histogram approach
nBinSize <- 5000 # number of bins
ggplot() +
    geom_histogram(data = subset(methylation_slim, strand == "+"),
                   aes(x = start, y = after_stat(count), fill = strand),
                   alpha = 0.5, bins = nBinSize) +
    geom_histogram(data = subset(methylation_slim, strand == "-"),
                   aes(x = start, y = -after_stat(count), fill = strand),
                   alpha = 0.5, bins = nBinSize) +
    scale_fill_manual(values = c("+" = "blue", "-" = "red")) +
    facet_grid(treatment ~ .) +
    labs(
        title = paste0("Distribution of Methylation Sites by Strand ",clone_of_interest, "_category", COI),
        x = "Chromosome Position",
        y = "Count (+ strand up, - strand down)",
        fill = "Strand"
    ) +
    theme_minimal() +
    geom_hline(yintercept = 0, linetype = "dashed")

(fN <- paste("Methylation_hist_on_Chromosomev2_", clone_of_interest, "_category", COI, ".png", sep = ""))
# save the plot
ggsave(fN, width = 20, height = 6)


# Check what treatments we have
table(methylation_slim$treatment)

# Basic counts by treatment
methylation_slim %>%
    count(treatment, name = "total_sites") %>%
    arrange(desc(total_sites))

# Visualize total sites per treatment
methylation_slim %>%
    count(treatment) %>%
    ggplot(aes(x = treatment, y = n, fill = treatment)) +
    geom_col() +
    theme_minimal() +
    labs(title = "Total Methylation Sites by Treatment",
         x = "Treatment", y = "Number of Sites") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Strand preferences by treatment
methylation_slim %>%
    count(treatment, strand) %>%
    ggplot(aes(x = treatment, y = n, fill = strand)) +
    geom_col(position = "dodge") +
    theme_minimal() +
    labs(title = "Strand Distribution by Treatment",
         x = "Treatment", y = "Number of Sites")
#save plot
(fN <- paste("Methylation_strand_by_Treatment", clone_of_interest, "_category", COI, ".png", sep = ""))
ggsave(fN, width = 10, height = 6)

# Which loci are affected differently by treatments?
locus_treatment_summary <- methylation_slim %>%
    count(locus_tag, treatment) %>%
    pivot_wider(names_from = treatment, values_from = n, values_fill = 0)

# Show loci with biggest differences between treatments
head(locus_treatment_summary) # take out here regions with no counts
# Filter for loci with significant differences

# loci not methylated
notMethylatedInPASN <- locus_treatment_summary %>%
    filter(PASN == 0) |> select(locus_tag) |> distinct()

notMethylatedInTSB <- locus_treatment_summary %>%
    filter(TSB == 0) |> select(locus_tag) |> distinct()

# join in AGU from desc
# bring in AGU
desc <- readRDS("descFasta.rds")
# write desc out
#write_tsv(desc, file = "all_GenesdescFasta.tsv")
#

# chromosome visualization
# Prepare data for chromosome visualization
methylation_plot_data <- methylation_slim %>%
    mutate(
        # Create y-position based on strand and treatment
        y_pos = case_when(
            strand == "+" & treatment == unique(methylation_slim$treatment)[1] ~ 1,
            strand == "-" & treatment == unique(methylation_slim$treatment)[1] ~ -1,
            strand == "+" & treatment == unique(methylation_slim$treatment)[2] ~ 0.5,
            strand == "-" & treatment == unique(methylation_slim$treatment)[2] ~ -0.5,
            TRUE ~ 0
        ),
        # Create treatment-strand combination for coloring
        treatment_strand = paste(treatment, strand, sep = "_")
    )

# Create the chromosome plot
chromosome_plot <- ggplot(methylation_plot_data, aes(x = start, y = y_pos)) +
    # Add horizontal reference lines for each treatment
    geom_hline(yintercept = c(1, 0.5, -0.5, -1),
               color = "lightgray", linetype = "dashed", alpha = 0.5) +

    # Add methylation sites as vertical ticks
    geom_segment(aes(xend = start, yend = 0, color = treatment_strand),
                 linewidth = 0.5, alpha = 0.7) +

    # Customize colors for treatment-strand combinations
    scale_color_manual(
        values = setNames(
            c("#E31A1C", "#FF7F00", "#1F78B4", "#33A02C"),
            c(paste(unique(methylation_slim$treatment)[1], "+", sep = "_"),
              paste(unique(methylation_slim$treatment)[1], "-", sep = "_"),
              paste(unique(methylation_slim$treatment)[2], "+", sep = "_"),
              paste(unique(methylation_slim$treatment)[2], "-", sep = "_"))
        ),
        name = "Treatment & Strand"
    ) +

    # Customize the plot
    scale_y_continuous(
        breaks = c(-1, -0.5, 0.5, 1),
        labels = c(paste(unique(methylation_slim$treatment)[1], "(-)"),
                   paste(unique(methylation_slim$treatment)[2], "(-)"),
                   paste(unique(methylation_slim$treatment)[2], "(+)"),
                   paste(unique(methylation_slim$treatment)[1], "(+)"))
    ) +

    theme_minimal() +
    theme(
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 10),
        legend.position = "bottom"
    ) +

    labs(
        title = "Methylation Sites Across Chromosome",
        subtitle = "Vertical ticks show methylation sites, + strand up, - strand down",
        x = "Genomic Position (bp)",
        y = "Treatment & Strand"
    )

# Display the plot
print(chromosome_plot)

#Optional: Create a zoomed-in version for a specific region
#Uncomment and adjust the range as needed
zoom_plot <- chromosome_plot +
    coord_cartesian(xlim = c(1.2E6, 1.4E6)) +
    labs(title = "Methylation Sites - Zoomed View")

print(zoom_plot)
# save plot
(fN <- paste("zoomed_Methylation_Chromosome_Plots_", clone_of_interest, "_category", COI, ".png", sep = ""))
ggsave(fN, width = 10, height = 6)

# Summary statistics
cat("Summary of methylation sites by treatment and strand:\n")
methylation_slim %>%
    count(treatment, strand) %>%
    arrange(treatment, strand) %>%
    print()

# Density comparison
cat("\nMethylation density (sites per kb) by treatment:\n")
methylation_slim %>%
    group_by(treatment) %>%
    summarise(
        min_pos = min(start),
        max_pos = max(start),
        total_sites = n(),
        genome_span_kb = (max_pos - min_pos) / 1000,
        density_per_kb = total_sites / genome_span_kb
    ) %>%
    print()


# visualize more

# Get treatment names
treatments <- unique(methylation_slim$treatment)
print(paste("Treatments found:", paste(treatments, collapse = ", ")))

# 1. SIDE-BY-SIDE CHROMOSOME PLOTS (FACETED)
# This separates treatments into different panels for easy comparison
faceted_plot <- methylation_slim %>%
    mutate(
        y_pos = ifelse(strand == "+", 1, -1),
        strand_label = paste("Strand", strand)
    ) %>%
    ggplot(aes(x = start, y = y_pos, color = strand)) +
    geom_segment(aes(xend = start, yend = 0),
                 size = 0.6, alpha = 0.8) +
    geom_hline(yintercept = 0, color = "black", size = 0.3) +
    facet_grid(treatment ~ ., scales = "free_y") +
    scale_color_manual(values = c("+" = "#E31A1C", "-" = "#1F78B4")) +
    scale_y_continuous(breaks = c(-1, 1), labels = c("- strand", "+ strand")) +
    theme_minimal() +
    theme(
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        strip.text = element_text(size = 12, face = "bold")
    ) +
    labs(
        title = "Methylation Sites by Treatment - Side by Side Comparison",
        x = "Genomic Position (bp)",
        y = "Strand",
        color = "Strand"
    )

print(faceted_plot)

zoom_plot2 <- faceted_plot +
    coord_cartesian(xlim = c(1.2E6, 1.4E6)) +
    labs(title = "Methylation Sites - Zoomed View")

print(zoom_plot2)
# save plot
(fN <- paste("Methylation_SideBySide_Chromosome_Plots_", clone_of_interest, "_category", COI, ".png", sep = ""))
ggsave(fN, width = 10, height = 6)

# 2. DENSITY PLOTS - Shows methylation density along the chromosome
density_plot <- methylation_slim %>%
    ggplot(aes(x = start, fill = treatment, color = treatment)) +
    geom_density(alpha = 0.6, size = 0.1, adjust = smoothin_param) +
    facet_wrap(~ strand, ncol = 1, labeller = label_both) +
    scale_fill_manual(values = c("#E31A1C", "#1F78B4")) +
    scale_color_manual(values = c("#E31A1C", "#1F78B4")) +
    theme_minimal() +
    labs(
        title = "Methylation Density Distribution by Treatment and Strand",
        x = "Genomic Position (bp)",
        y = "Density",
        fill = "Treatment",
        color = "Treatment"
    )

print(density_plot)
(fN <- paste("Methylation_DensityOverlay_on_Chromosomev2_", clone_of_interest, "_category", COI, ".png", sep = ""))
# save the plot
ggsave(fN, width = 20, height = 6)



# 3. DIFFERENTIAL METHYLATION - Shows regions where treatments differ
# Create genomic bins to compare treatments
bin_size <- 5000  # 10kb bins
methylation_binned <- methylation_slim %>%
    mutate(bin = floor(start / bin_size) * bin_size) %>%
    count(treatment, strand, bin, name = "sites_in_bin") %>%
    complete(treatment, strand, bin, fill = list(sites_in_bin = 0)) %>%
    pivot_wider(names_from = treatment, values_from = sites_in_bin, values_fill = 0)

# find the peak region
methylation_binned
methylation_binned$TSBnPASN <- methylation_binned$TSB + methylation_binned$PASN
methylation_binned$endBin <- methylation_binned$bin + bin_size
methylation_binned[order(methylation_binned$TSBnPASN, decreasing = TRUE)[1:10],]

# bin_size <- 5000
# looking for all features! # check this directly on Biocyc Genome Viewer
# A tibble: 10 × 6
# strand     bin  PASN   TSB TSBnPASN  endBin
# <chr>    <dbl> <int> <int>    <int>   <dbl>
#     1 -      1865000   196   168      364 1870000
# 2 +       505000    79    71      150  510000
# 3 +      1255000    66    64      130 1260000
# 4 -      2145000    65    62      127 2150000
# 5 +       810000    57    51      108  815000
# 6 -      2025000    46    49       95 2030000
# 7 +      1865000    48    47       95 1870000
# 8 -      1910000    42    46       88 1915000
# 9 +       465000    42    39       81  470000
# 10 -      1255000    39    41       80 1260000


# Calculate differences between treatments
# Get the actual treatment names from the original data
treatment_names <- unique(methylation_slim$treatment)
treatment1 <- treatment_names[1]
treatment2 <- treatment_names[2]

# Check if we have exactly 2 treatments
if(length(treatment_names) != 2) {
    stop(paste("Expected 2 treatments, found:", length(treatment_names)))
}

methylation_binned <- methylation_binned %>%
    mutate(
        difference = !!sym(treatment1) - !!sym(treatment2),
        abs_difference = abs(difference),
        dominant_treatment = case_when(
            difference > 0 ~ treatment1,
            difference < 0 ~ treatment2,
            TRUE ~ "Equal"
        )
    )

# Plot the differences
diff_plot <- methylation_binned %>%
    filter(abs_difference > 0) %>%  # Only show bins with differences
    ggplot(aes(x = bin, y = difference, fill = dominant_treatment)) +
    geom_col(alpha = 0.8) +
    geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
    facet_wrap(~ strand, ncol = 1, labeller = label_both) +
    scale_fill_manual(values = c("#E31A1C", "#1F78B4", "gray50")) +
    theme_minimal() +
    labs(
        title = paste("Differential Methylation Between", treatment1, "vs", treatment2),
        subtitle = paste("Positive values =", treatment1, "higher; Negative values =", treatment2, "higher"),
        x = "Genomic Position (bp)",
        y = paste("Difference in Sites (", treatment1, "-", treatment2, ")"),
        fill = "Higher in"
    )

print(diff_plot)
# save plot
(fN <- paste("Methylation_Differential_Plots_", clone_of_interest, "_category", COI, ".png", sep = ""))
ggsave(fN, width = 20, height = 6)

# figure out the one max loci
methylation_binned[order(methylation_binned$abs_difference, decreasing = TRUE),]

# one peak with 10 more in pasn -> where is it
maxBin_PASN <- as.numeric(methylation_binned[order(methylation_binned$abs_difference, decreasing = TRUE)[1],2])


# what genes are in this one
# prepare desc
desc$strand <- "+"
desc$strand[grepl("complement", desc$location)] <- "-"
table(desc$strand)
desc$location_1 <- gsub(x = gsub(x = desc$location, pattern = "complement\\(", replacement = ""), pattern = "\\)", replacement = "")
desc$geneStart <- as.numeric(sapply(strsplit(desc$location_1, split = "\\.\\."), function(x)x[1]))
desc$geneEnd <- as.numeric(sapply(strsplit(desc$location_1, split = "\\.\\."), function(x)x[2]))
desc$location_1 <- NULL

write_tsv(desc, file = "all_GenesdescFasta.tsv")
#
maxBinStart <- maxBin_PASN
maxBinEnd <- maxBin_PASN+bin_size

# Filter rows where any part of the gene is within the region
subset_proteins <- desc[(desc$geneStart <= maxBinEnd & desc$geneEnd >= maxBinStart),]
print(subset_proteins)

# empty here? -> yes no protein coding genes in this region.
# really no genes?
# biocyc: https://biocyc.org/tu?orgid=GCF_000462955&id=TU2BZV-949
# https://biocyc.org/genbro/genbro.shtml?orgid=GCF_000462955&replicon=NC_022222
# region is full of t-RNAs

# now we want to look at the TSB region
# another interesting one 6 more in TSB 8 vs 2
maxBin_TSB <- as.numeric(methylation_binned[order(methylation_binned$abs_difference, decreasing = TRUE)[2],2])

maxBin_TSB
maxBinStart <- maxBin_TSB
maxBinEnd <- maxBin_TSB+bin_size

# Filter rows where any part of the gene is within the region
subset_proteins <- desc[(desc$geneStart <= maxBinEnd & desc$geneEnd >= maxBinStart),]
print(subset_proteins)

#
# another interesting one 14 vs 19 -> 5 more in TSB
maxBin_TSB <- as.numeric(methylation_binned[order(methylation_binned$abs_difference, decreasing = TRUE)[4],2])

maxBin_TSB
maxBinStart <- maxBin_TSB
maxBinEnd <- maxBin_TSB+bin_size

# Filter rows where any part of the gene is within the region
subset_proteins <- desc[(desc$geneStart <= maxBinEnd & desc$geneEnd >= maxBinStart),]
print(subset_proteins)
# translocase and sporulation protein


# another interesting one 18 vs 23 -> 5 more in TSB
maxBin_TSB <- as.numeric(methylation_binned[order(methylation_binned$abs_difference, decreasing = TRUE)[6],2])

maxBin_TSB
maxBinStart <- maxBin_TSB
maxBinEnd <- maxBin_TSB+bin_size

# Filter rows where any part of the gene is within the region
subset_proteins <- desc[(desc$geneStart <= maxBinEnd & desc$geneEnd >= maxBinStart),]
print(subset_proteins)
# battery of ribosomal proteins
# what protein expression do these show?

plot(prX$diff.PASNvsTSB_givenAncestor[prX$AGU %in% subset_proteins$proteinID], type="l")

# how many differences of methylations do we see in the binned df
methylation_binned %>%
    filter(abs_difference > 0) %>%
    summarise(total_differences = n())



# Hunt the motifs with more fine granular bins
# 3. DIFFERENTIAL METHYLATION - Shows regions where treatments differ
# Create genomic bins to compare treatments
bin_size <- 1  # 2 bp bins
methylation_fine_binned <- methylation_slim %>%
    mutate(bin = floor(start / bin_size) * bin_size) %>%
    count(treatment, strand, bin, name = "sites_in_bin") %>%
    complete(treatment, strand, bin, fill = list(sites_in_bin = 0)) %>%
    pivot_wider(names_from = treatment, values_from = sites_in_bin, values_fill = 0)

# Calculate differences between treatments
# Get the actual treatment names from the original data
treatment_names <- unique(methylation_slim$treatment)
treatment1 <- treatment_names[1]
treatment2 <- treatment_names[2]

# Check if we have exactly 2 treatments
if(length(treatment_names) != 2) {
    stop(paste("Expected 2 treatments, found:", length(treatment_names)))
}

methylation_fine_binned <- methylation_fine_binned %>%
    mutate(
        difference = !!sym(treatment1) - !!sym(treatment2),
        abs_difference = abs(difference),
        dominant_treatment = case_when(
            difference > 0 ~ treatment1,
            difference < 0 ~ treatment2,
            TRUE ~ "Equal"
        )
    )

# Plot the differences
# diff_fine_plot <- methylation_fine_binned %>%
#     filter(abs_difference > 0) %>%  # Only show bins with differences
#     ggplot(aes(x = bin, y = difference, fill = dominant_treatment)) +
#     geom_col(alpha = 0.8) +
#     geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
#     facet_wrap(~ strand, ncol = 1, labeller = label_both) +
#     scale_fill_manual(values = c("#E31A1C", "#1F78B4", "gray50")) +
#     theme_minimal() +
#     labs(
#         title = paste("Differential Methylation Between", treatment1, "vs", treatment2),
#         subtitle = paste("Positive values =", treatment1, "higher; Negative values =", treatment2, "higher"),
#         x = "Genomic Position (bp)",
#         y = paste("Difference in Sites (", treatment1, "-", treatment2, ")"),
#         fill = "Higher in"
#     )
# print(diff_fine_plot)
# # save plot
# (fN <- paste("Methylation_Differential_Fine_Plots_", clone_of_interest, "_category", COI, ".png", sep = ""))
# ggsave(fN, width = 20, height = 6)

# how many differences of methylations do we see in the binned df
methylation_fine_binned %>%
    filter(abs_difference > 0) %>%
    summarise(total_differences = n())


# join with MtX_wPrX to find Motif
MtX_wPrX_slim <- mtX_wPrX %>%
    select(locus_tag, start, end, strand, feature, treatment, group, category,
           base_before, base_after, motif, gene_start, gene_end, AGU, mean_iBAQ,
           parsedProteinDescription, diff.PASNvsTSB_givenAncestor,
           FDR.PASNvsTSB_givenAncestor, statistic.PASNvsTSB_givenAncestor) %>%
    distinct()

# join in the binned data
methylation_differences <- methylation_fine_binned %>%
    filter(abs_difference > 0)

methylation_fine_binned_wMtX <- left_join(methylation_differences, MtX_wPrX_slim, by = c("bin" = "start", "strand" = "strand", "dominant_treatment" = "treatment"))

plot(methylation_fine_binned_wMtX$difference, methylation_fine_binned_wMtX$diff.PASNvsTSB_givenAncestor)
boxplot(methylation_fine_binned_wMtX$diff.PASNvsTSB_givenAncestor ~ methylation_fine_binned_wMtX$difference )


# go for logos
library(ggseqlogo)

# Assuming your dataframe is named `df`

# Motifs with positive difference
motifs_pos <- methylation_fine_binned_wMtX$motif[methylation_fine_binned_wMtX$difference > 0]

# Motifs with negative difference
motifs_neg <- methylation_fine_binned_wMtX$motif[methylation_fine_binned_wMtX$difference < 0]

# Function to convert a character vector of sequences into a character matrix
motif_to_matrix <- function(motifs) {
    seqs <- strsplit(motifs, split = "")
    do.call(rbind, seqs)
}

mat_pos <- motif_to_matrix(motifs_pos)
mat_neg <- motif_to_matrix(motifs_neg)


# Just pass the vector of sequences
ggseqlogo(motifs_neg, method = 'prob') + ggtitle("Negative Difference Motifs")
ggseqlogo(motifs_pos, method = 'prob') + ggtitle("Positive Difference Motifs")



PASN_meth <- methylation_fine_binned_wMtX |> filter(difference > 0) |> select(AGU) |> group_by(AGU) |> summarise(uniqueMethylations = n())
TSB_meth <- methylation_fine_binned_wMtX |> filter(difference < 0) |> select(AGU) |> group_by(AGU) |> summarise(uniqueMethylations = n())

# find exclusive AGUs in lists
PASN_meth <- PASN_meth$AGU
TSB_meth <- TSB_meth$AGU
# find exclusive AGUs in lists
exclusive_PASN <- setdiff(PASN_meth, TSB_meth)
exclusive_TSB <- setdiff(TSB_meth, PASN_meth)

# write out the data
(fN <- paste("PASNspecific_methylation_", clone_of_interest, "_category", COI, ".tsv", sep = ""))
write.table(exclusive_PASN, file = fN,quote = FALSE, row.names = FALSE) # boring
# GO:0005622
# Intracellular anatomical structure
# 180 of 993	0.12	0.42	0.00047
# GO:0005737
# Cytoplasm
# 166 of 905	0.12	0.42	0.00059
# GO:0110165
# Cellular anatomical entity
# 313 of 2035	0.04	0.4	0.00019


(fN <- paste("TSBspecific_methylation_", clone_of_interest, "_category", COI, ".tsv", sep = ""))
write.table(exclusive_TSB, file = fN,quote = FALSE, row.names = FALSE) # boring


# go for distance from start and end AND perc. of total
methylation_fine_binned_wMtX$distance_from_start <- methylation_fine_binned_wMtX$bin - methylation_fine_binned_wMtX$gene_start
methylation_fine_binned_wMtX$distance_from_end <- methylation_fine_binned_wMtX$gene_end - methylation_fine_binned_wMtX$bin
# calculate percentage of total gene length
methylation_fine_binned_wMtX$gene_length <- methylation_fine_binned_wMtX$gene_end - methylation_fine_binned_wMtX$gene_start
methylation_fine_binned_wMtX$perc_from_start <- methylation_fine_binned_wMtX$distance_from_start / methylation_fine_binned_wMtX$gene_length * 100

# visualize
# ggplot(methylation_slim, aes(x = start)) +
#     geom_density(data = subset(methylation_slim, strand == "+"),
#                  aes(y = after_stat(density), fill = strand),
#                  alpha = 0.5, adjust = smoothin_param) +  # <-- less smoothing) +
#     geom_density(data = subset(methylation_slim, strand == "-"),
#                  aes(y = -after_stat(density), fill = strand),
#                  alpha = 0.5,
#                  adjust = smoothin_param) +  # <-- less smoothing) +
#     scale_fill_manual(values = c("+" = "blue", "-" = "red")) +
#     facet_grid(treatment ~ .) +
#     labs(
#         title = "Density of Methylation Sites by Strand",
#         x = "Chromosome Position",
#         y = "Density",
#         fill = "Strand"
#     ) +
#     theme_minimal() +
#     geom_hline(yintercept = 0, linetype = "dashed")

# use density and absolute distance from start
ggplot(methylation_fine_binned_wMtX, aes(x = distance_from_start)) +
    geom_density(aes(fill = dominant_treatment), alpha = 0.5, adjust = smoothin_param) +
    scale_fill_manual(values = c("#E31A1C", "#1F78B4", "gray50")) +
    labs(
        title = "Density of Methylation Sites by Distance from Start",
        x = "Distance from Start (bp)",
        y = "Density",
        fill = "Dominant Treatment"
    ) +
    theme_minimal()


# use density and absolute distance from end
ggplot(methylation_fine_binned_wMtX, aes(x = distance_from_end)) +
    geom_density(aes(fill = dominant_treatment), alpha = 0.5, adjust = smoothin_param) +
    scale_fill_manual(values = c("#E31A1C", "#1F78B4", "gray80")) +
    labs(
        title = "Density of Methylation Sites by Distance from End",
        x = "Distance from End (bp)",
        y = "Density",
        fill = "Dominant Treatment"
    ) +
    theme_minimal()


# visualize percentage from start but separate for treatment PASN up and TSB down
ggplot(methylation_fine_binned_wMtX, aes(x = perc_from_start)) +
    geom_density(aes(fill = dominant_treatment), alpha = 0.5, adjust = smoothin_param) +
    scale_fill_manual(values = c("#E31A1C", "#1F78B4", "gray50")) +
    labs(
        title = "Density of Methylation Sites by Percentage from Start",
        x = "Percentage from Start (%)",
        y = "Density",
        fill = "Dominant Treatment"
    ) +
    theme_minimal()


# any useful logo for motifs within the first 50bp
# Motifs with positive difference
motifs_pos_inFirstX <- methylation_fine_binned_wMtX$motif[methylation_fine_binned_wMtX$difference > 0 & methylation_fine_binned_wMtX$distance_from_start <= 20]

# Motifs with negative difference
motifs_neg_inFirstX <- methylation_fine_binned_wMtX$motif[methylation_fine_binned_wMtX$difference < 0 & methylation_fine_binned_wMtX$distance_from_start <= 20]

# Just pass the vector of sequences
pdf("motifs_inFirstX.pdf", width = 10, height = 10)
ggseqlogo(motifs_neg, method = 'prob') + ggtitle("Negative Difference Motifs")
ggseqlogo(motifs_pos, method = 'prob') + ggtitle("Positive Difference Motifs")
ggseqlogo(motifs_neg_inFirstX, method = 'prob') + ggtitle("Negative Difference Motifs in First X")
ggseqlogo(motifs_pos_inFirstX, method = 'prob') + ggtitle("Positive Difference Motifs in First X")
dev.off()


