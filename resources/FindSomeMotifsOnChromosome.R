#!/usr/local/bin/Rscript
#
# Jonas Grossmann <jg@fgcz.ethz.ch>
# 2024
#
#

# here we would like to find motifs on the full chromosome

# GATC Motif Density Analysis for Staphylococcus aureus Chromosome
# Load required libraries
library(ggplot2)
library(dplyr)
library(seqinr)
library(tidyr)
# go for logos
library(ggseqlogo)


# Your chromosome sequence (remove header and combine lines)
# Note: This is a simplified example sequence for demonstration purposes.
# sequence <- "ATTAACTTGTGGATAATTATTAACATGGTGTGTTTAGAAGTTATCCACGGCTGTTATTTTTGTGTATAACTTAAAAATTTAAGAAAGATGGAGTAAATTTATGTCGGAAAAAGAAATTTGGGAAAAAGTGCTTGAAATTGCTCAAGAAAAATTATCAGCTGTAAGTTACTCAACTTTCCTAAAAGATACTGAGCTTTACACGATCAAAGATGGTGAAGCTATCGTATTATCGAGTATTCCTTTTAATGCAAATTGGTTAAATCAACAATATGCTGAAATTATCCAAGCAATCTTATTTGATGTTGTAGGCTATGAAGTAAAACCTCACTTTATTACTACTGAAGAATTAGCAAATTATAGTAATAATGAAACTGCTACTCCAAAAGAAACAACAAAACCTTCTACTGAAACAACTGAGGATAATCATGTGCTTGGTAGAGAGCAATTCAATGCCCATAACACATTTGACACTTTTGTAATCGGACCTGGTAACCGCTTTCCACATGCAGCGAGTTTAGCTGTGGCCGAAGCACCAGCCAAAGCGTACAATCCATTATTTTATCTATGGAGGTGTTGGTTTAGGAAAAACCCATTTAATGCATGCCATTGGTCATCATGTTTTAGATAATAA"

# read sequence from a file
my_SA_DNA <- read.fasta(file = "SA_6850_ncbi_sequence.fasta",
                        seqtype = "DNA",
                        as.string = TRUE,
                        forceDNAtolower = FALSE)

str(my_SA_DNA)  # Check the structure of the loaded sequences)

sequence <- my_SA_DNA$CP006706.1  # Assuming the first sequence in the file is the chromosome

# Function to find all GATC motif positions
find_gatc_positions <- function(seq, my_motif = "GATC") {
    positions <- c()
    start <- 1

    while (TRUE) {
        pos <- regexpr(my_motif, substr(seq, start, nchar(seq)))
        if (pos == -1) break

        actual_pos <- start + pos - 1
        positions <- c(positions, actual_pos)
        start <- actual_pos + 1
    }

    return(positions)
}

# Find GATC positions
gatc_positions <- find_gatc_positions(sequence)
cat("Found", length(gatc_positions), "GATC motifs\n")
cat("Positions:", gatc_positions, "\n")

# Create a data frame for plotting
gatc_df <- data.frame(position = gatc_positions)

# Calculate sequence length
seq_length <- nchar(sequence)
cat("Sequence length:", seq_length, "bp\n")

# Create density plot
p1 <- ggplot(gatc_df, aes(x = position)) +
    geom_density(fill = "blue", alpha = 0.3, color = "darkblue") +
    xlim(0, seq_length) +
    labs(title = "GATC Motif Density Distribution",
         subtitle = paste("Staphylococcus aureus chromosome (", length(gatc_positions), "GATC sites)"),
         x = "Chromosome Position (bp)",
         y = "Density") +
    theme_minimal()

print(p1)

# Alternative: Create histogram with bins (useful for overlaying with methylation data)
p2 <- ggplot(gatc_df, aes(x = position)) +
    geom_histogram(bins = 50, fill = "blue", alpha = 0.7, color = "darkblue") +
    xlim(0, seq_length) +
    labs(title = "GATC Motif Distribution (Histogram)",
         subtitle = paste("Staphylococcus aureus chromosome (", length(gatc_positions), "GATC sites)"),
         x = "Chromosome Position (bp)",
         y = "Count") +
    theme_minimal()

print(p2)

smoothin_param <- 1E-6
p4 <- ggplot(gatc_df, aes(x = position)) +
    geom_density(aes(y = after_stat(density)), fill = "blue", alpha = 0.5, adjust = smoothin_param)
print(p4)




# Create binned data that you can use for overlaying with methylation data
# Adjust bin_size as needed to match your methylation data bins
bin_size <- 1000
seq_length <- nchar(sequence)
n_bins <- ceiling(seq_length / bin_size)

# Create binned data
# q I want to generate binned data from my gatc positions
gatc_binned <- data.frame(
    bin = 1:(n_bins-1),
    gatc_count = integer(n_bins-1),
    bin_start = seq(1, (seq_length-bin_size), by = bin_size),
    bin_end = seq(bin_size, seq_length, by = bin_size)
)

# Calculate GATC counts per bin
for (i in 1:(n_bins-1)) {
    bin_start <- gatc_binned$bin_start[i]
    bin_end <- gatc_binned$bin_end[i]
    gatc_binned$gatc_count[i] <- sum(gatc_positions >= bin_start & gatc_positions < bin_end)
}


gatc_binned$bin_center <- (gatc_binned$bin_start + gatc_binned$bin_end) / 2

# Show the binned data structure
head(gatc_binned, 10)

# Plot binned GATC data
p3 <- ggplot(gatc_binned, aes(x = bin_center, y = gatc_count)) +
    geom_col(fill = "blue", alpha = 0.7, width = bin_size * 0.8) +
    labs(title = paste("GATC Motifs per", bin_size, "bp Bin"),
         subtitle = "Ready for overlay with methylation data",
         x = "Chromosome Position (bp)",
         y = "GATC Count per Bin") +
    theme_minimal()

print(p3)


smoothin_param <- 1E-6
p4 <- ggplot(gatc_binned, aes(x = bin_center, y = gatc_count)) +
    geom_density(aes(y = after_stat(density)), fill = "blue", alpha = 0.5, adjust = smoothin_param)

print(p4)





# do we have motifs anyway in our motifs?
mtX <- readRDS("methylation_data_EGGnogAnnotated.rds")

## we want to filter out modified_base -> it is unclear what type of feature it is
mtX <- mtX[mtX$feature != "modified_base",]
# anyway we are only looking into ancestor (the rest is for students of Bio253)
mtX <- mtX[mtX$group == "6850",]

table(mtX$feature)

mtX_A <- mtX[mtX$feature == "m6A",]
mtX_C <- mtX[mtX$feature == "m4C",]

#



# Motifs with positive difference
motifs_A <- mtX_A |> select(start, feature, motif) |> distinct() |> pull(motif)

# Motifs with negative difference
# motifs_C <- mtX_C$motif
motifs_C <- mtX_C |> select(start, feature, motif) |> distinct() |> pull(motif)

# Just pass the vector of sequences
pdf("All_Motifs_6mA_4mC.pdf", width = 10, height = 5)
ggseqlogo(motifs_A, method = 'prob') + ggtitle("Motifs 6mA")
ggseqlogo(motifs_C, method = 'prob') + ggtitle("Motifs 4mC")
dev.off()


# for particular category
# Motifs with positive difference
whatCategory <- "gene"
motifs_A <- mtX_A |> select(start, feature, motif, category) |> filter(category == whatCategory) |> distinct() |> pull(motif)
motifs_C <- mtX_C |> select(start, feature, motif, category) |> filter(category == whatCategory) |> distinct() |> pull(motif)


# Just pass the vector of sequences
pdf(paste0("All_Motifs_Category_",whatCategory,".pdf"), width = 20, height = 5)
ggseqlogo(motifs_A, method = 'bits') + ggtitle(paste("Motifs for Category: ", whatCategory, " 6mA")) + theme_logo()
ggseqlogo(motifs_C, method = 'bits') + ggtitle(paste("Motifs for Category: ", whatCategory, " 4mC")) + theme_logo()
dev.off()



# Just pass the vector of sequences
pdf(paste0("All_Motifs_Category_",whatCategory,".pdf"), width = 20, height = 5)
ggseqlogo(motifs_A, method = 'prob', col_scheme='base_pairing') + ggtitle(paste("Motifs for Category: ", whatCategory, " 6mA")) + theme_logo()
ggseqlogo(motifs_C, method = 'prob', col_scheme='base_pairing') + ggtitle(paste("Motifs for Category: ", whatCategory, " 4mC")) + theme_logo()
dev.off()

