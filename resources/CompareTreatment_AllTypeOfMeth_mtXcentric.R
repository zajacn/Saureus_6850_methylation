# 2025-05-28
# this script is used to join the prX data with the methylation data
# storing the results in the resources folder is too big

# load libraries
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)

# read in the data
prX <- readRDS("SA6850_prXallWide_moreMetaInfo.rds")
mtX <- readRDS("methylation_data_EGGnogAnnotated.rds")

# here mtX centric
# join PrX 2 mtx
mtX_wPrX <- left_join(x = mtX, y = prX, by = c("locus_tag" = "LocusTag"))




# what catetories are there
table(mtX_wPrX$category)

# category of Interest
COI <- "downstream"


# we want to filter out modified_base -> it is unclear what type of feature it is
mtX_wPrX <- mtX_wPrX[mtX_wPrX$feature != "modified_base",]
table(mtX_wPrX$feature)
table(mtX_wPrX$category)

# filter for Ancestor only
unique(mtX_wPrX$group)

# anyway we are only looking into ancestor
mtX_wPrX <- mtX_wPrX[mtX_wPrX$group == "6850",]

# we are only interested in one type features
mtX_wPrX <- mtX_wPrX[mtX_wPrX$category == COI,]
table(mtX_wPrX$feature)
table(mtX_wPrX$category)

# we are only interested in the TSB treatment
condition_of_interest <- "TSB"
mtX_wPrX_OneTreatment <- mtX_wPrX[mtX_wPrX$treatment == condition_of_interest,]

# now summarize per AGU and treatment
table(mtX_wPrX_OneTreatment$feature)
table(mtX_wPrX_OneTreatment$category)




# Hi Jonas,
#
# Thanks a lot. I am not sure if using log2-ratio for the methylation data is the way to go as a single methylation could already affect transcription. How I would move forward:
#     1.       Focus on Ancestor_TSB only
# a.       Describe the global methylation pattern of SA6850 (what is methylated, how is it methylated, hyper vs hypo methylation,6mA vs 4mC)
# 2. Eliminate all genes were we have the exact same methylations (Anc_TSB vs Anc_SN) -> Open question: Are methylation always happening at the exact same position within a gene? Or could the same amount of methylation but on different positions already influence regulations?
#     3. Have a list methylated regions that are different from Anc_TSB vs Anc_SN -> do we even find this?
#     4. If yes, match with proteomics and see if we find regions of different methylation count and significant proteomics.
#
# Hope this makes sense.
# Best, Lukas

#     1.       Focus on Ancestor_TSB only
# a.       Describe the global methylation pattern of SA6850 (what is methylated, how is it methylated, hyper vs hypo methylation,6mA vs 4mC)

# slim
methylation_slim <- mtX_wPrX %>% select(treatment, start, strand, feature, locus_tag) |> distinct() # here we do have multiple reps in

# for the plotting
clone_of_interest <- "Ancestor 6850"



# visualize data as density on the chromosome with start as position
# Create the density plot with proper direction based on strand
ggplot(methylation_slim, aes(x = start)) +
    geom_density(data = subset(methylation_slim, strand == "+"),
                 aes(y = after_stat(density), fill = strand),
                 alpha = 0.5, adjust = 0.01) +  # <-- less smoothing) +
    geom_density(data = subset(methylation_slim, strand == "-"),
                 aes(y = -after_stat(density), fill = strand),
                 alpha = 0.5,
                 adjust = 0.01) +  # <-- less smoothing) +
    scale_fill_manual(values = c("+" = "blue", "-" = "red")) +
    labs(
        title = "Density of m6A Methylation Sites by Strand",
        x = "Chromosome Position",
        y = "Density",
        fill = "Strand"
    ) +
    theme_minimal() +
    geom_hline(yintercept = 0, linetype = "dashed")

# stat density
ggplot(methylation_slim, aes(x = start)) +
    stat_density(data = subset(methylation_slim, strand == "+"),
                 aes(y = after_stat(density), fill = strand),
                 alpha = 0.5) +
    stat_density(data = subset(methylation_slim, strand == "-"),
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
                   alpha = 0.5, bins = 1500) +
    geom_histogram(data = subset(methylation_slim, strand == "-"),
                   aes(x = start, y = -after_stat(count), fill = strand),
                   alpha = 0.5, bins = 1500) +
    scale_fill_manual(values = c("+" = "blue", "-" = "red")) +
    labs(
        title = paste0("Distribution of Methylation Sites by Strand (",clone_of_interest," ",condition_of_interest,")"),
        x = "Chromosome Position",
        y = "Count (+ strand up, - strand down)",
        fill = "Strand"
    ) +
    theme_minimal() +
    geom_hline(yintercept = 0, linetype = "dashed")

(fN <- paste("Methylation_hist_on_Chromosome_allFeatures", clone_of_interest, "_", condition_of_interest, ".png", sep = ""))
# save the plot
ggsave(fN, width = 10, height = 6)


# now we want to summarize the methylation data
str(methylation_slim)

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

# Which loci are affected differently by treatments?
locus_treatment_summary <- methylation_slim %>%
    count(locus_tag, treatment) %>%
    pivot_wider(names_from = treatment, values_from = n, values_fill = 0)

# Show loci with biggest differences between treatments
head(locus_treatment_summary) # take out here regions with no counts
# Filter for loci with significant differences

#
notMethylatedInPASN <- locus_treatment_summary %>%
    filter(PASN == 0) |> select(locus_tag) |> distinct()

notMethylatedInTSB <- locus_treatment_summary %>%
    filter(TSB == 0) |> select(locus_tag) |> distinct()

# join in AGU from desc
# bring in AGU
desc <- readRDS("descFasta.rds")
# write desc out
write_tsv(desc, file = "all_GenesdescFasta.tsv")

notMethylatedInPASN <- left_join(notMethylatedInPASN, desc, by = c("locus_tag" = "LocusTag"))
notMethylatedInTSB <- left_join(notMethylatedInTSB, desc, by = c("locus_tag" = "LocusTag"))




# write to file
write_tsv(x = notMethylatedInPASN |> select(proteinID), file = "notMethylatedInPASN_AGU.tsv")
write_tsv(x = notMethylatedInTSB |> select(proteinID), file = "notMethylatedInTSB_AGU.tsv")

write_tsv(x = notMethylatedInPASN, file = "notMethylatedInPASN.tsv")
write_tsv(x = notMethylatedInTSB, file = "notMethylatedInTSB.tsv")


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
        subtitle = "Vertical ticks show m4C sites, + strand up, - strand down",
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
library(dplyr)
library(ggplot2)
library(tidyr)

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

# 2. DENSITY PLOTS - Shows methylation density along the chromosome
density_plot <- methylation_slim %>%
    ggplot(aes(x = start, fill = treatment, color = treatment)) +
    geom_density(alpha = 0.6, size = 1) +
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

# 3. DIFFERENTIAL METHYLATION - Shows regions where treatments differ
# Create genomic bins to compare treatments
bin_size <- 10000  # 10kb bins
methylation_binned <- methylation_slim %>%
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

# TREATMENT-SPECIFIC SITES
# Find sites that are unique or highly enriched in one treatment
treatment_names <- unique(methylation_slim$treatment)
treatment_comparison <- methylation_slim %>%
    count(start, strand, treatment) %>%
    pivot_wider(names_from = treatment, values_from = n, values_fill = 0) %>%
    mutate(
        total_sites = rowSums(select(., all_of(treatment_names))),
        treatment1_specific = !!sym(treatment_names[1]) > 0 & !!sym(treatment_names[2]) == 0,
        treatment2_specific = !!sym(treatment_names[1]) == 0 & !!sym(treatment_names[2]) > 0,
        shared = !!sym(treatment_names[1]) > 0 & !!sym(treatment_names[2]) > 0,
        site_type = case_when(
            treatment1_specific ~ paste(treatment_names[1], "specific"),
            treatment2_specific ~ paste(treatment_names[2], "specific"),
            shared ~ "Shared",
            TRUE ~ "None"
        )
    )

# Plot treatment-specific sites
specific_plot <- treatment_comparison %>%
    filter(site_type != "None") %>%
    ggplot(aes(x = start, y = strand, color = site_type)) +
    geom_point(alpha = 0.7, size = 1.5) +
    scale_color_manual(values = c("#E31A1C", "#1F78B4", "#FF7F00")) +
    theme_minimal() +
    labs(
        title = "Treatment-Specific vs Shared Methylation Sites",
        x = "Genomic Position (bp)",
        y = "Strand",
        color = "Site Type"
    )

print(specific_plot)

# 5. SUMMARY STATISTICS
cat("\n=== TREATMENT COMPARISON SUMMARY ===\n")

# Overall counts
cat("\nTotal sites per treatment:\n")
methylation_slim %>%
    count(treatment, sort = TRUE) %>%
    print()

# Strand distribution
cat("\nStrand distribution by treatment:\n")
methylation_slim %>%
    count(treatment, strand) %>%
    pivot_wider(names_from = strand, values_from = n, values_fill = 0) %>%
    mutate(strand_ratio = `+` / `-`) %>%
    print()

# Treatment-specific summary
cat("\nTreatment-specific sites:\n")
treatment_comparison %>%
    count(site_type) %>%
    mutate(percentage = round(n / sum(n) * 100, 1)) %>%
    print()

# Genomic coverage
cat("\nGenomic span by treatment:\n")
methylation_slim %>%
    group_by(treatment) %>%
    summarise(
        min_pos = min(start),
        max_pos = max(start),
        span_kb = (max_pos - min_pos) / 1000,
        sites = n(),
        density_per_kb = sites / span_kb
    ) %>%
    print()




# Now let's look at GAPs that are different between treatments.
library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)

# Function to calculate gaps between methylation sites
calculate_gaps <- function(data, treatment_name, strand_name) {
    sites <- data %>%
        filter(treatment == treatment_name, strand == strand_name) %>%
        arrange(start) %>%
        pull(start)

    if(length(sites) <= 1) {
        return(data.frame(
            treatment = treatment_name,
            strand = strand_name,
            gap_start = numeric(0),
            gap_end = numeric(0),
            gap_length = numeric(0)
        ))
    }

    gaps <- data.frame(
        treatment = treatment_name,
        strand = strand_name,
        gap_start = sites[-length(sites)],
        gap_end = sites[-1],
        gap_length = sites[-1] - sites[-length(sites)]
    )

    return(gaps)
}

# Calculate gaps for all treatment-strand combinations
treatments <- unique(methylation_slim$treatment)
strands <- unique(methylation_slim$strand)

all_gaps <- map_dfr(treatments, function(t) {
    map_dfr(strands, function(s) {
        calculate_gaps(methylation_slim, t, s)
    })
})

print(paste("Total gaps found:", nrow(all_gaps)))
print("Gap summary by treatment:")
all_gaps %>% count(treatment, sort = TRUE)

# 1. GAP LENGTH DISTRIBUTION COMPARISON
gap_distribution <- all_gaps %>%
    ggplot(aes(x = gap_length, fill = treatment)) +
    geom_histogram(alpha = 0.7, bins = 50, position = "identity") +
    facet_wrap(~ strand, scales = "free_y", labeller = label_both) +
    scale_x_log10(labels = scales::comma) +
    scale_fill_manual(values = c("#E31A1C", "#1F78B4")) +
    theme_minimal() +
    labs(
        title = "Distribution of Gap Lengths Between Treatments",
        x = "Gap Length (bp, log scale)",
        y = "Count",
        fill = "Treatment"
    )

print(gap_distribution)

# 2. GAP LENGTHS AS CHROMOSOME VIEW
# Show gaps as blocks along the chromosome
gap_chromosome <- all_gaps %>%
    filter(gap_length >= quantile(gap_length, 0.25)) %>%  # Focus on larger gaps
    mutate(
        y_pos = case_when(
            strand == "+" & treatment == treatments[1] ~ 1,
            strand == "-" & treatment == treatments[1] ~ -1,
            strand == "+" & treatment == treatments[2] ~ 0.5,
            strand == "-" & treatment == treatments[2] ~ -0.5
        )
    ) %>%
    ggplot(aes(y = y_pos, color = treatment)) +
    geom_segment(aes(x = gap_start, xend = gap_end, yend = y_pos),
                 size = 2, alpha = 0.7) +
    geom_hline(yintercept = c(-1, -0.5, 0.5, 1),
               color = "lightgray", linetype = "dashed", alpha = 0.5) +
    scale_color_manual(values = c("#E31A1C", "#1F78B4")) +
    scale_y_continuous(
        breaks = c(-1, -0.5, 0.5, 1),
        labels = c(paste(treatments[1], "(-)"),
                   paste(treatments[2], "(-)"),
                   paste(treatments[2], "(+)"),
                   paste(treatments[1], "(+)"))
    ) +
    theme_minimal() +
    theme(
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank()
    ) +
    labs(
        title = "Large Unmethylated Gaps Along Chromosome",
        subtitle = "Thick lines show gaps without methylation sites",
        x = "Genomic Position (bp)",
        y = "Treatment & Strand",
        color = "Treatment"
    )

print(gap_chromosome)

# 3. TREATMENT-SPECIFIC GAPS
# Find gaps that exist in one treatment but not the other
gap_comparison <- all_gaps %>%
    # Create overlapping windows to find similar gaps
    mutate(
        gap_center = (gap_start + gap_end) / 2,
        gap_window = round(gap_center / 5000) * 5000  # 5kb windows
    ) %>%
    group_by(strand, gap_window) %>%
    summarise(
        treatments_present = list(unique(treatment)),
        n_treatments = length(unique(treatment)),
        avg_gap_length = mean(gap_length),
        min_start = min(gap_start),
        max_end = max(gap_end),
        .groups = "drop"
    ) %>%
    mutate(
        gap_type = case_when(
            n_treatments == 1 & treatments[1] %in% unlist(treatments_present) ~ paste(treatments[1], "specific"),
            n_treatments == 1 & treatments[2] %in% unlist(treatments_present) ~ paste(treatments[2], "specific"),
            n_treatments == 2 ~ "Shared",
            TRUE ~ "Other"
        )
    )

# Plot treatment-specific gaps
specific_gaps_plot <- gap_comparison %>%
    filter(gap_type != "Other") %>%
    ggplot(aes(x = (min_start + max_end)/2, y = avg_gap_length, color = gap_type)) +
    geom_point(alpha = 0.7, size = 2) +
    facet_wrap(~ strand, labeller = label_both) +
    scale_color_manual(values = c("#E31A1C", "#1F78B4", "#FF7F00")) +
    scale_y_log10(labels = scales::comma) +
    theme_minimal() +
    labs(
        title = "Treatment-Specific vs Shared Unmethylated Gaps",
        x = "Genomic Position (bp)",
        y = "Average Gap Length (bp, log scale)",
        color = "Gap Type"
    )

print(specific_gaps_plot)

# 4. LARGE GAP ANALYSIS
# Focus on the biggest gaps and see how they differ
large_gaps_threshold <- quantile(all_gaps$gap_length, 0.9)  # Top 10% of gaps

large_gaps <- all_gaps %>%
    filter(gap_length >= large_gaps_threshold) %>%
    arrange(desc(gap_length))

large_gaps_plot <- large_gaps %>%
    ggplot(aes(x = gap_start, y = gap_length, color = treatment, shape = strand)) +
    geom_point(size = 3, alpha = 0.8) +
    scale_color_manual(values = c("#E31A1C", "#1F78B4")) +
    scale_y_log10(labels = scales::comma) +
    theme_minimal() +
    labs(
        title = paste("Largest Unmethylated Gaps (>", scales::comma(large_gaps_threshold), "bp)"),
        x = "Gap Start Position (bp)",
        y = "Gap Length (bp, log scale)",
        color = "Treatment",
        shape = "Strand"
    )

print(large_gaps_plot)

# 5. STATISTICAL SUMMARY
cat("\n=== GAP ANALYSIS SUMMARY ===\n")

# Overall gap statistics
cat("\nGap length statistics by treatment:\n")
all_gaps %>%
    group_by(treatment) %>%
    summarise(
        n_gaps = n(),
        mean_gap = round(mean(gap_length)),
        median_gap = round(median(gap_length)),
        max_gap = max(gap_length),
        total_gap_length = sum(gap_length)
    ) %>%
    print()

# Strand-specific gap statistics
cat("\nGap statistics by treatment and strand:\n")
all_gaps %>%
    group_by(treatment, strand) %>%
    summarise(
        n_gaps = n(),
        mean_gap = round(mean(gap_length)),
        median_gap = round(median(gap_length)),
        .groups = "drop"
    ) %>%
    print()

# Treatment-specific gap regions
cat("\nTreatment-specific gap regions:\n")
gap_comparison %>%
    count(gap_type) %>%
    mutate(percentage = round(n / sum(n) * 100, 1)) %>%
    print()

# Largest gaps by treatment
cat("\nTop 10 largest gaps:\n")
large_gaps %>%
    select(treatment, strand, gap_start, gap_end, gap_length) %>%
    head(10) %>%
    print()

# invoke here new file to get locus_tags that are in these gaps.


# Statistical test for gap length differences
if(length(treatments) == 2) {
    gap1 <- all_gaps %>% filter(treatment == treatments[1]) %>% pull(gap_length)
    gap2 <- all_gaps %>% filter(treatment == treatments[2]) %>% pull(gap_length)

    wilcox_result <- wilcox.test(gap1, gap2)
    cat(sprintf("\nWilcoxon test for gap length differences between %s and %s:\n",
                treatments[1], treatments[2]))
    cat(sprintf("p-value: %.2e\n", wilcox_result$p.value))
    cat(sprintf("Interpretation: %s\n",
                ifelse(wilcox_result$p.value < 0.05, "Significant difference", "No significant difference")))
}



# invoke locus tag in these gaps
head(desc)
# significance = case_when(
#     FDR < 0.05 & diff > 1 ~ "Up",
#     FDR < 0.05 & diff < -1 ~ "Down",
#     TRUE ~ "Not Significant"
# )

desc$strand <- "+"
desc$strand[grepl("complement", desc$location)] <- "-"
table(desc$strand)
desc$location_1 <- gsub(x = gsub(x = desc$location, pattern = "complement\\(", replacement = ""), pattern = "\\)", replacement = "")
desc$geneStart <- as.numeric(sapply(strsplit(desc$location_1, split = "\\.\\."), function(x)x[1]))
desc$geneEnd <- as.numeric(sapply(strsplit(desc$location_1, split = "\\.\\."), function(x)x[2]))
desc$location_1 <- NULL
# now we want to join the gaps with the desc
#
# library(dplyr)
# library(ggplot2)
# library(tidyr)
# library(purrr)

# Function to find genes overlapping with gaps
find_overlapping_genes <- function(gaps_df, genes_df) {
    overlaps <- map_dfr(1:nrow(gaps_df), function(i) {
        gap <- gaps_df[i, ]

        # Find genes that overlap with this gap
        overlapping <- genes_df %>%
            filter(
                # Gene start or end is within the gap
                (geneStart >= gap$gap_start & geneStart <= gap$gap_end) |
                    (geneEnd >= gap$gap_start & geneEnd <= gap$gap_end) |
                    # Gap is entirely within the gene
                    (geneStart <= gap$gap_start & geneEnd >= gap$gap_end) |
                    # Gene spans the entire gap
                    (geneStart <= gap$gap_start & geneEnd >= gap$gap_end)
            )

        if(nrow(overlapping) > 0) {
            overlapping %>%
                mutate(
                    gap_start = gap$gap_start,
                    gap_end = gap$gap_end,
                    gap_length = gap$gap_length,
                    gap_treatment = gap$treatment,
                    gap_strand = gap$gap_strand,
                    overlap_type = case_when(
                        geneStart >= gap$gap_start & geneEnd <= gap$gap_end ~ "Gene within gap",
                        geneStart >= gap$gap_start & geneStart <= gap$gap_end ~ "Gene start in gap",
                        geneEnd >= gap$gap_start & geneEnd <= gap$gap_end ~ "Gene end in gap",
                        geneStart <= gap$gap_start & geneEnd >= gap$gap_end ~ "Gap within gene",
                        TRUE ~ "Other overlap"
                    )
                )
        } else {
            NULL
        }
    })

    return(overlaps)
}

# Make sure we have the gaps data (assuming it was calculated in previous code)
# If not, recalculate gaps
if(!exists("all_gaps")) {
    # Recalculate gaps function
    calculate_gaps <- function(data, treatment_name, strand_name) {
        sites <- data %>%
            filter(treatment == treatment_name, strand == strand_name) %>%
            arrange(start) %>%
            pull(start)

        if(length(sites) <= 1) {
            return(data.frame(
                treatment = treatment_name,
                gap_strand = strand_name,
                gap_start = numeric(0),
                gap_end = numeric(0),
                gap_length = numeric(0)
            ))
        }

        gaps <- data.frame(
            treatment = treatment_name,
            gap_strand = strand_name,
            gap_start = sites[-length(sites)],
            gap_end = sites[-1],
            gap_length = sites[-1] - sites[-length(sites)]
        )

        return(gaps)
    }

    treatments <- unique(methylation_slim$treatment)
    strands <- unique(methylation_slim$strand)

    all_gaps <- map_dfr(treatments, function(t) {
        map_dfr(strands, function(s) {
            calculate_gaps(methylation_slim, t, s)
        })
    })
}

# Rename strand column if needed
if("strand" %in% names(all_gaps)) {
    all_gaps <- all_gaps %>% rename(gap_strand = strand)
}

# Focus on only the TOP 10 LARGEST gaps per treatment
top_gaps <- all_gaps %>%
    group_by(treatment) %>%
    arrange(desc(gap_length)) %>%
    slice_head(n = 10) %>%
    ungroup()

# write out top_gaps
write_tsv(top_gaps, file = "top_10_largest_methylation_gaps_perTreatment.tsv")

cat("=== FOCUSING ON TOP 10 LARGEST GAPS PER TREATMENT ===\n")
cat("Top gaps by treatment:\n")
top_gaps %>%
    group_by(treatment) %>%
    summarise(
        n_gaps = n(),
        min_gap_size = min(gap_length),
        max_gap_size = max(gap_length),
        median_gap_size = median(gap_length)
    ) %>%
    print()

# Find genes overlapping with ONLY the top gaps
genes_in_gaps <- find_overlapping_genes(top_gaps, desc)
# write out
write_tsv(genes_in_gaps, file = "genes_in_top_methylation_gaps.tsv") # verified by jonas the first 2 gaps are correct still some caps are in PASN AND TSB


cat("=== GENES IN TOP METHYLATION GAPS ANALYSIS ===\n")
cat(paste("Top gaps analyzed:", nrow(top_gaps), "\n"))
cat(paste("Genes found in top gaps:", nrow(genes_in_gaps), "\n"))
cat(paste("Unique genes affected:", length(unique(genes_in_gaps$LocusTag)), "\n"))

# 1. SUMMARY OF GENES AFFECTED BY TREATMENT
genes_by_treatment <- genes_in_gaps %>%
    group_by(gap_treatment, LocusTag) %>%
    summarise(
        n_gaps = n(),
        total_gap_length = sum(gap_length),
        gene_description = first(description),
        gene_start = first(geneStart),
        gene_end = first(geneEnd),
        gene_length = gene_end - gene_start,
        .groups = "drop"
    )

cat("\nGenes affected by gaps in each treatment:\n")
genes_by_treatment %>%
    count(gap_treatment, name = "genes_affected") %>%
    print()

# some bu here
# 2. TREATMENT-SPECIFIC GENE EFFECTS
# Find genes that are in gaps in one treatment but not the other
treatment_specific_genes <- genes_by_treatment %>%
    select(LocusTag, gap_treatment, gene_description) %>%
    pivot_wider(names_from = gap_treatment, values_from = gap_treatment,
                values_fill = "absent") %>%
    mutate(
        gene_gap_pattern = case_when(
            get(names(.)[2]) != "absent" & get(names(.)[3]) == "absent" ~ paste("Only", names(.)[2]),
            get(names(.)[2]) == "absent" & get(names(.)[3]) != "absent" ~ paste("Only", names(.)[3]),
            get(names(.)[2]) != "absent" & get(names(.)[3]) != "absent" ~ "Both treatments",
            TRUE ~ "Neither"
        )
    )

cat("\nTreatment-specific gene gap patterns:\n")
treatment_specific_genes %>%
    count(gene_gap_pattern) %>%
    print()

# there is a but in this table
treatment_specific_genes
# one column does not make sense
treatment_specific_genes$gene_gap_pattern <- NULL


# 3. VISUALIZATION: TOP GAPS AND THEIR GENES
top_gaps_plot <- top_gaps %>%
    ggplot(aes(x = gap_start, y = gap_treatment, color = gap_treatment)) +
    geom_segment(aes(xend = gap_end, yend = gap_treatment), size = 4, alpha = 0.8) +
    geom_text(aes(x = (gap_start + gap_end)/2, label = paste0(round(gap_length/1000, 1), "kb")),
              vjust = -0.5, size = 3) +
    scale_color_manual(values = c("#E31A1C", "#1F78B4")) +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 12)) +
    labs(
        title = "Top 10 Largest Unmethylated Gaps per Treatment",
        x = "Genomic Position (bp)",
        y = "Treatment",
        color = "Treatment"
    )

print(top_gaps_plot)

# Plot genes within these top gaps
if(nrow(genes_in_gaps) > 0) {
    genes_gaps_plot <- genes_in_gaps %>%
        ggplot(aes(x = geneStart, y = gap_treatment, color = gap_treatment)) +
        geom_segment(aes(xend = geneEnd, yend = gap_treatment), size = 3, alpha = 0.7) +
        geom_point(aes(x = (geneStart + geneEnd)/2), size = 2) +
        scale_color_manual(values = c("#E31A1C", "#1F78B4")) +
        theme_minimal() +
        theme(axis.text.y = element_text(size = 12)) +
        labs(
            title = "Genes Located Within Top 10 Largest Gaps",
            x = "Genomic Position (bp)",
            y = "Treatment",
            color = "Treatment"
        )

    print(genes_gaps_plot)
}

# 4. GAP SIZE VS GENE SIZE ANALYSIS
gap_gene_sizes <- genes_in_gaps %>%
    mutate(
        gene_length = geneEnd - geneStart,
        gap_to_gene_ratio = gap_length / gene_length
    )

size_comparison_plot <- gap_gene_sizes %>%
    ggplot(aes(x = gene_length, y = gap_length, color = gap_treatment)) +
    geom_point(alpha = 0.7, size = 2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
    scale_x_log10(labels = scales::comma) +
    scale_y_log10(labels = scales::comma) +
    scale_color_manual(values = c("#E31A1C", "#1F78B4")) +
    theme_minimal() +
    labs(
        title = "Gap Length vs Gene Length",
        subtitle = "Dashed line shows equal gap and gene lengths",
        x = "Gene Length (bp, log scale)",
        y = "Gap Length (bp, log scale)",
        color = "Treatment"
    )

print(size_comparison_plot)

# 5. FUNCTIONAL ANALYSIS OF AFFECTED GENES
# Look for patterns in gene descriptions
cat("\n=== FUNCTIONAL ANALYSIS ===\n")

# Most common gene types in gaps
cat("\nMost frequently affected gene types (by description keywords):\n")
gene_keywords <- genes_in_gaps %>%
    select(LocusTag, description, gap_treatment) %>%
    distinct() %>%
    mutate(
        description_lower = tolower(description),
        gene_type = case_when(
            grepl("ribosom|rRNA|tRNA", description_lower) ~ "RNA/Ribosomal",
            grepl("transport|permease|ABC", description_lower) ~ "Transport",
            grepl("transcription|regulat|repressor", description_lower) ~ "Regulation",
            grepl("DNA|replication|repair", description_lower) ~ "DNA processes",
            grepl("metaboli|enzyme|synthase|dehydrogenase", description_lower) ~ "Metabolism",
            grepl("membrane|cell wall|peptidoglycan", description_lower) ~ "Cell structure",
            grepl("hypothetical|unknown|conserved", description_lower) ~ "Hypothetical/Unknown",
            TRUE ~ "Other"
        )
    )

gene_keywords %>%
    count(gene_type, gap_treatment) %>%
    pivot_wider(names_from = gap_treatment, values_from = n, values_fill = 0) %>%
    arrange(desc(rowSums(select(., -gene_type)))) %>%
    print()

# 6. DETAILED GENE LIST
cat("\n=== TOP AFFECTED GENES ===\n")
cat("Genes with largest gaps (top 20):\n")
genes_in_gaps %>%
    arrange(desc(gap_length)) %>%
    select(LocusTag, description, gap_treatment, gap_strand, gap_length, geneStart, geneEnd) %>%
    head(20) %>%
    print()

# 7. TREATMENT-SPECIFIC AFFECTED GENES
treatments <- unique(genes_in_gaps$gap_treatment)
if(length(treatments) == 2) {
    cat(sprintf("\nGenes affected ONLY in %s:\n", treatments[1]))
    treatment_specific_genes %>%
        filter(gene_gap_pattern == paste("Only", treatments[1])) %>%
        select(LocusTag, gene_description) %>%
        head(10) %>%
        print()

    cat(sprintf("\nGenes affected ONLY in %s:\n", treatments[2]))
    treatment_specific_genes %>%
        filter(gene_gap_pattern == paste("Only", treatments[2])) %>%
        select(LocusTag, gene_description) %>%
        head(10) %>%
        print()
}

# 9. DETAILED LIST OF TOP GAPS WITH GENES
cat("\n=== DETAILED TOP GAPS WITH GENES ===\n")
top_gaps_with_genes <- top_gaps %>%
    left_join(
        genes_in_gaps %>%
            select(gap_start, gap_end, gap_treatment, LocusTag, description, geneStart, geneEnd, overlap_type),
        by = c("gap_start", "gap_end", "treatment" = "gap_treatment")
    ) %>%
    arrange(treatment, desc(gap_length))

for(t in unique(top_gaps$treatment)) {
    cat(sprintf("\n--- TOP GAPS IN %s ---\n", t))

    gaps_this_treatment <- top_gaps_with_genes %>%
        filter(treatment == t) %>%
        arrange(desc(gap_length))

    for(i in 1:nrow(gaps_this_treatment)) {
        gap <- gaps_this_treatment[i, ]
        cat(sprintf("\nGap #%d: %s - %s (%s bp)\n",
                    i,
                    scales::comma(gap$gap_start),
                    scales::comma(gap$gap_end),
                    scales::comma(gap$gap_length)))

        if(!is.na(gap$LocusTag)) {
            cat(sprintf("  Gene: %s - %s\n", gap$LocusTag, gap$description))
            cat(sprintf("  Overlap: %s\n", gap$overlap_type))
        } else {
            cat("  No genes found in this gap\n")
        }
    }
}
# 10. FINAL SUMMARY
cat("\n=== FINAL SUMMARY OF GAP ANALYSIS ===\n")
cat(sprintf("Total gaps analyzed: %d\n", nrow(all_gaps)))
cat(sprintf("Total unique genes affected: %d\n", length(unique(genes_in_gaps$LocusTag))))
cat(sprintf("Total gaps with genes: %d\n", nrow(genes_in_gaps)))
cat(sprintf("Total gaps in top 10 largest gaps: %d\n", nrow(top_gaps)))
cat(sprintf("Total genes in top gaps: %d\n", nrow(genes_in_gaps)))



# integrate PrX data for Ancestor
str(prX)

prX_anc <- prX %>%
    select(protein_Id, nrPeptides, IDcolumn, LocusTag, protein_length,
           mean_iBAQ, AGU, gene, diff.PASNvsTSB_givenAncestor,
           FDR.PASNvsTSB_givenAncestor, statistic.PASNvsTSB_givenAncestor)



# working on Excel file to get the locus tags that are in the gaps but not in other condition
# Success in StringDB
# Looking at LocusTags in Gaps that are mutually exclusive for treatment
# GAP in PASN condition and NOT in TSB
# enrichment for polysaccharide
# https://version-12-0.string-db.org/cgi/network?networkId=bFP7Q9tZ2TFs
# textmining evidence: https://version-12-0.string-db.org/cgi/textmining?networkId=bKDVllAwVR0Y

# GAPs in TSB condition and NOT in PASN
# enrichement
# https://version-12-0.string-db.org/cgi/network?networkId=bGbyFXDO9pCy



PASNspecificGenesInGaps <- treatment_specific_genes |> filter(PASN=="PASN", TSB=="absent") |> select(LocusTag, gene_description)
TSBspecificGenesInGaps <- treatment_specific_genes |> filter(PASN=="absent", TSB=="TSB") |> select(LocusTag, gene_description)


# join with prX_anc
PASNspecific_proteinExpression <- PASNspecificGenesInGaps |> left_join(prX_anc, by = "LocusTag")
TSBspecific_proteinExpression <- TSBspecificGenesInGaps |> left_join(prX_anc, by = "LocusTag")
# write out
write_tsv(PASNspecific_proteinExpression, file = "PASN_specificGaps_withProteinExpressionData.tsv")
write_tsv(TSBspecific_proteinExpression, file = "TSB_specificGaps_withProteinExpressionData.tsv")

# Visualize the Diffs as barplots
library(ggplot2)
# Barplot for PASN-specific genes in gaps
p1 <- ggplot(PASNspecific_proteinExpression, aes(x = reorder(LocusTag, diff.PASNvsTSB_givenAncestor), y = diff.PASNvsTSB_givenAncestor)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    labs(title = "PASN-Specific Genes in Gaps",
         x = "Locus Tag",
         y = "Difference (PASN - TSB)") +
    theme_minimal()
# Barplot for TSB-specific genes in gaps
p2 <- ggplot(TSBspecific_proteinExpression, aes(x = reorder(LocusTag, diff.PASNvsTSB_givenAncestor), y = diff.PASNvsTSB_givenAncestor)) +
    geom_bar(stat = "identity", fill = "darkorange") +
    coord_flip() +
    labs(title = "TSB-Specific Genes in Gaps",
         x = "Locus Tag",
         y = "Difference (PASN - TSB)") +
    theme_minimal()
# Print the plots
print(p1)
print(p2)

# generate input files for StringDB
stringTSB_gappedProteins <- TSBspecific_proteinExpression %>%
    select(AGU, diff.PASNvsTSB_givenAncestor) |> filter(!is.na(AGU))
write_tsv(stringTSB_gappedProteins, file = "TSB_specific_gapped_proteins.tsv", col_names = FALSE)

stringPASN_gappedProteins <- PASNspecific_proteinExpression %>%
    select(AGU, diff.PASNvsTSB_givenAncestor) |> filter(!is.na(AGU))
write_tsv(stringPASN_gappedProteins, file = "PASN_specific_gapped_proteins.tsv", col_names = FALSE)

hist(PASNspecific_proteinExpression$diff.PASNvsTSB_givenAncestor, breaks = 50, main = "PASN-Specific Gap Differences", xlab = "Difference (PASN - TSB)", col = "steelblue")
hist(TSBspecific_proteinExpression$diff.PASNvsTSB_givenAncestor, breaks = 50, main = "TSB-Specific Gap Differences", xlab = "Difference (PASN - TSB)", col = "darkorange")

# stringNW
# PASN specific genes in gaps
# https://version-12-0.string-db.org/cgi/network?networkId=b7IpYyVSFwHn

# TSB specific genes in gaps
# https://version-12-0.string-db.org/cgi/network?networkId=bO2HUdvD4PMN



