#!/usr/local/bin/Rscript
#
# Jonas Grossmann <jg@fgcz.ethz.ch>
# 2025
#
#
# 2025-06-25
# this script is used to join the prX data with the methylation data

# load libraries
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(purrr)
library(rtracklayer)

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
COI <- "allCat_6mAonly_10bp" # this will also be the description of the files outputed


# we want to filter out modified_base -> it is unclear what type of feature it is
mtX_wPrX <- mtX_wPrX[mtX_wPrX$feature != "modified_base",]
table(mtX_wPrX$feature)
mtX_wPrX <- mtX_wPrX[mtX_wPrX$feature != "m4C",]
#mtX_wPrX <- mtX_wPrX[mtX_wPrX$feature != "m6A",]
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


# make df more manageable -> still in here there might be doubs
mtX_wPrX %>% select(treatment, start, strand, feature, locus_tag) |> distinct() |> dim()
# slim
methylation_slim <- mtX_wPrX %>% select(treatment, start, strand, feature, locus_tag) |> distinct() # here we do have multiple reps in

# for the plotting
clone_of_interest <- "Ancestor6850"


smoothin_param <- 0.001

# do it individually for m6A and m4C
# visualize data as density on the chromosome with start as position
# Create the density plot with proper direction based on strand
ggplot(methylation_slim, aes(x = start)) +
    geom_density(data = subset(methylation_slim, strand == "+"),
                 aes(y = after_stat(density), fill = feature),
                 alpha = 0.5, adjust = smoothin_param) +  # <-- less smoothing) +
    geom_density(data = subset(methylation_slim, strand == "-"),
                 aes(y = -after_stat(density), fill = feature),
                 alpha = 0.5,
                 adjust = smoothin_param) +  # <-- less smoothing) +
    scale_fill_manual(values = c("m6A" = "blue", "m4C" = "yellow")) +
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

# Check what treatments we have
table(methylation_slim$treatment)


# 3. DIFFERENTIAL METHYLATION - Shows regions where treatments differ
# Create genomic bins to compare treatments
bin_size <- 10  # 10bp bins
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
    geom_col(alpha = 0.99) +
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
ggsave(fN, width = 40, height = 6)

# write to table
(fn <- paste("Methylation_binned_", clone_of_interest, "_category_", COI, ".tsv", sep = ""))
write_tsv(methylation_binned, fn)


# loot into the gff to find features in bins?
mygff <- readGFF("SA_6850_GCA_000462955.1_ASM46295v1_genomic.gff", version = 3)
colnames(mygff)

mygffdf <- as.data.frame(mygff)

# 2025-06-18
# Find genomic features that overlap with bins having abs_difference > 2
# First, filter methylation_binned for high difference bins
high_diff_bins <- methylation_binned %>%
    filter(abs_difference > 2)

# Function to check if a genomic feature overlaps with a bin
# A feature overlaps if its start OR end position falls within the bin range
find_overlapping_features <- function(methylation_bins, gff_data) {

    overlapping_features <- data.frame()

    for (i in 1:nrow(methylation_bins)) {
        bin_start <- methylation_bins$bin[i]
        bin_end <- methylation_bins$endBin[i]
        bin_strand <- methylation_bins$strand[i]

        # Convert strand notation: "-" becomes "-", anything else becomes "+"
        gff_strand_match <- if(bin_strand == "-") "-" else "+"

        # Find features where:
        # 1. The strand matches
        # 2. The feature's start OR end position falls within the bin range
        matches <- gff_data %>%
            filter(
                strand == gff_strand_match,
                (start >= bin_start & start <= bin_end) |
                    (end >= bin_start & end <= bin_end)
            ) %>%
            mutate(
                methylation_bin_start = bin_start,
                methylation_bin_end = bin_end,
                methylation_strand = bin_strand,
                PASN_count = methylation_bins$PASN[i],
                TSB_count = methylation_bins$TSB[i],
                difference = methylation_bins$difference[i],
                abs_difference = methylation_bins$abs_difference[i],
                dominant_treatment = methylation_bins$dominant_treatment[i]
            )

        overlapping_features <- rbind(overlapping_features, matches)
    }

    return(overlapping_features)
}

# Alternative approach using dplyr joins for better performance with large datasets
library(dplyr)

# Create expanded bin ranges for joining
expanded_bins <- high_diff_bins %>%
    rowwise() %>%
    do(data.frame(
        bin_pos = seq(.$bin, .$endBin, by = 1),
        methylation_strand = .$strand,
        PASN_count = .$PASN,
        TSB_count = .$TSB,
        difference = .$difference,
        abs_difference = .$abs_difference,
        dominant_treatment = .$dominant_treatment,
        bin_start = .$bin,
        bin_end = .$endBin
    )) %>%
    ungroup()

# More efficient approach: direct overlap detection
efficient_overlap <- function(methylation_bins, gff_data) {
    results <- list()

    for (i in 1:nrow(methylation_bins)) {
        bin_start <- methylation_bins$bin[i]
        bin_end <- methylation_bins$endBin[i]
        bin_strand <- methylation_bins$strand[i]

        # Convert strand
        gff_strand_match <- if(bin_strand == "-") "-" else "+"

        # Find overlapping features
        overlaps <- gff_data %>%
            filter(
                strand == gff_strand_match,
                # Check if feature overlaps with bin (feature start <= bin end AND feature end >= bin start)
                start <= bin_end & end >= bin_start
            ) %>%
            mutate(
                bin_start = bin_start,
                bin_end = bin_end,
                methylation_strand = bin_strand,
                PASN_count = methylation_bins$PASN[i],
                TSB_count = methylation_bins$TSB[i],
                difference = methylation_bins$difference[i],
                abs_difference = methylation_bins$abs_difference[i],
                dominant_treatment = methylation_bins$dominant_treatment[i]
            )

        if (nrow(overlaps) > 0) {
            results[[i]] <- overlaps
        }
    }

    if (length(results) > 0) {
        return(do.call(rbind, results))
    } else {
        return(data.frame())
    }
}


# do it for all methylation marks
all_results <- efficient_overlap(methylation_binned, mygffdf)


# joining in PrX
all_results$TSBnPASN <- all_results$TSB_count + all_results$PASN_count

# Join with prX data
all_results_wPrX <- left_join(all_results, prX, by = c("locus_tag" = "LocusTag"))
(fn <- paste("All_Methylation_inBins_wPrX_", clone_of_interest, "_category_", COI, ".tsv", sep = ""))
write_tsv(all_results_wPrX, fn)


# from a friend:
# 🧬 1. Functional Importance of rRNA and tRNA Genes
# These loci are essential and highly conserved, playing central roles in protein synthesis.
# Methylation in these regions is often regulatory, either protecting important sequences from cleavage or helping modulate transcription and translation efficiency.
#
# 🧪 2. Protection from Restriction Enzymes (Restriction-Modification Systems)
# Bacteria use methylation as part of their defense system to distinguish self from foreign DNA.
# tRNA and rRNA genes are targeted for protective methylation to prevent cleavage by their own restriction enzymes.
# Type I, II, and III restriction-modification systems often target sequence motifs that may be abundant in rRNA/tRNA regions due to their conserved structure.
#
# 🧫 3. Transcriptional Regulation
# Methylation in these loci can modulate transcription under different growth or stress conditions.
# For example, 6mA methylation in bacteria has been associated with:
#     Gene silencing or activation
# Transcriptional pausing or enhancement
# Temporal control of ribosome production, especially during nutrient limitation or stationary phase.
#
# 🧬 4. Sequence Bias and Methylation Motifs
# The motifs recognized by S. aureus DNA methyltransferases may coincidentally be overrepresented in tRNA/rRNA gene regions.
# You can check this by:
#     Identifying known methylation motifs (e.g., from PacBio SMRT or Nanopore data)
# Mapping their distribution across the genome
# Comparing it to GC content and functional regions
#
# 📈 5. High Transcriptional Activity Leads to Accessibility
# rRNA and tRNA genes are among the most actively transcribed in bacteria.
# Active transcription may leave DNA more accessible to methyltransferases, resulting in higher methylation density.
#
# ✅ Suggested Analysis
# To confirm and understand this better:
#
#     Map known methylation motifs across your genome and overlap them with tRNA/rRNA loci.
# Check if methylation is strand-specific or symmetric.
# Compare with gene expression data if available — are these methylation-rich loci also highly transcribed?
#     Compare different growth conditions — do methylation patterns in these regions change?
#
#


