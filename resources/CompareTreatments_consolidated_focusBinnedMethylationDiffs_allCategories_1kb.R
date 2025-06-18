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
COI <- "all_1kb" # this will also be the description of the files outputed


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


# make df more manageable -> still in here there might be doubs
mtX_wPrX %>% select(treatment, start, strand, feature, locus_tag) |> distinct() |> dim()
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
nBinSize <- 1000 # bin size
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

# Summary statistics
cat("Summary of methylation sites by treatment and strand:\n")
methylation_slim %>%
    count(treatment, strand) %>%
    arrange(treatment, strand) %>%
    print()



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
bin_size <- 1000  # 10kb bins
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

# write to table
(fn <- paste("Methylation_binned_", clone_of_interest, "_category_", COI, ".tsv", sep = ""))
write_tsv(methylation_binned, fn)
