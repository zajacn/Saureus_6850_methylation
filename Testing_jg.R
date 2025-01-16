#!/usr/local/bin/Rscript
#
# Jonas Grossmann <jg@fgcz.ethz.ch>
# 2024
#
#

# Jonas looks into methylation data
library(tidyverse)
library(ggplot2)
library(reshape2)

mtx <- readRDS("resources/methylation_data.rds")
dim(mtx)


# try to visualize mA6 as binary heatmap along the genome
str(mtx)
scoreFilter <- 2.5
strandFilter <- "+"

mtx_6mAplus <- mtx |> filter(feature == "m6A" & IPDRatio > scoreFilter & strand == strandFilter) |>
    select(c(Name, treatment, group, start)) |>
    mutate(value = 1)

# 2 conditions, 2 clones, ancestor
unique(mtx_6mAplus$Name)

dim(table(mtx_6mAplus$start, mtx_6mAplus$Name))
str(mtx_6mAplus)
mtx_6mAplus$condition <- paste(mtx_6mAplus$treatment, mtx_6mAplus$group, sep = "_")
# now plotting as a heatmap
p1 <- ggplot(mtx_6mAplus, aes(x = factor(start), y = factor(Name), fill = factor(value))) +
    geom_tile() +
    scale_fill_manual(values = c("white", "black")) +  # Black and White for binary
    #scale_fill_manual(values = c("0" = "white", "1" = "black")) +  # Ensure pure black and white
    labs(x = "", y = "Name", title = "Binary Heatmap: Names vs Start") +
    theme_minimal() +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, size=0.1),
        axis.text.y = element_text(size = 8)
    )

    # Color the columns based on Treatment_Group (e.g., differentiate columns)
    #facet_wrap(~condition, scales = "free_x", ncol = 1)  # Optional if you have different treatments
p1

# not satisfying
pdf("BinaryHeatmap.pdf", width = 40, height = 20)
print(p1)
dev.off()


p2 <- ggplot(mtx_6mAplus, aes(x = start, y = factor(Name), fill = factor(value))) +
    geom_tile() +
    #scale_fill_manual(values = c("white", "black")) +  # Black and White for binary
    scale_fill_manual(values = c("0" = "white", "1" = "black")) +  # Ensure pure black and white
    labs(x = "", y = "Name", title = "Binary Heatmap: Names vs Start") +
    theme_minimal() +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, size=0.1),
        axis.text.y = element_text(size = 8)
    )
p2

binary_data <- mtx_6mAplus %>%
    distinct(Name, start) %>%
    mutate(value = 1)

# Generate a matrix of all combinations of Name and Start
all_combinations <- expand.grid(
    Name = unique(binary_data$Name),
    start = unique(binary_data$start)
)

# max start value
mxStart <- max(all_combinations$start)

# Merge with the binary data to get 1's and 0's
final_data <- merge(all_combinations, binary_data, by = c("Name", "start"), all.x = TRUE)
final_data$value[is.na(final_data$value)] <- 0  # Replace NA with 0

# Plot the heatmap using ggplot2 with a continuous x-axis
p3 <- ggplot(final_data, aes(x = start, y = factor(Name), fill = factor(value))) +
    geom_tile() +
    #scale_fill_manual(values = c("0" = "white", "1" = "black")) +  # Black and White for binary
    scale_fill_manual(values = c("0" = "white", "1" = "black")) +  # Ensure pure black and white
    labs(x = "Start", y = "Name", title = "Binary Heatmap: Names vs Start") +
    theme_minimal() +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(size = 8),
        panel.grid = element_blank()  # Remove gridlines for cleaner look
    ) +
    scale_x_continuous(breaks = seq(min(final_data$start), max(final_data$start), by = 10))  # Customize x-axis breaks if needed

print(p3)



# use xtabs or table function
table(mtx_6mAplus$start, mtx_6mAplus$Name)
myXtab <- xtabs(~start + Name, data = mtx_6mAplus)
str(myXtab)

myXtab <- as.data.frame.matrix(myXtab)
myBoolmat <- myXtab > 0
dim(myXtab)
image(myBoolmat)

# replace TRUE with 1 and FALSE with 0 in my matrix
dim(myXtab) # here we also do have 2 in

myXtab[myXtab > 0] <- 1
myMat <- as.matrix(myXtab)
colnames(myMat)
rownames(myMat) <- as.numeric(rownames(myXtab))

# now plotting as a heatmap
myMat <- t(myMat)
dim(myMat)
heatmap(myMat, Rowv = NA, Colv = NA, col = c("white", "black"), scale = "none", margins = c(5, 5))

# Convert the matrix to long format (melt)

myMat_long <- melt(myMat)

# Rename columns for clarity
colnames(myMat_long) <- c("Name", "Start", "value")

# Plot the heatmap using ggplot2
ggplot(myMat_long, aes(x = Start, y = factor(Name), fill = factor(value))) +
    geom_tile() +
    scale_fill_manual(values = c("0" = "white", "1" = "black")) +  # Black and White for binary
    labs(x = "Start (Column Names)", y = "Name (Rows)", title = "Binary Heatmap: Names vs Start") +
    theme_minimal() +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1),  # Rotate x-axis labels for readability
        axis.text.y = element_text(size = 8),  # Adjust y-axis text size if needed
        panel.grid = element_blank()  # Remove gridlines for cleaner look
    ) +
    scale_x_discrete(breaks = colnames(myMat))  # Ensure x-axis corresponds to the column names



# from https://stackoverflow.com/questions/67227496/how-to-plot-timeseries-events-in-a-barcode-style

dim(myMat_long)
myMat_long$Name <- as.factor(myMat_long$Name)
myMat_long$Start <- as.numeric(myMat_long$Start)
str(myMat_long)

# Plot the heatmap using ggplot2
myChromosomePlot_plusStrand <- ggplot(myMat_long) + geom_point(aes(x = Start, y = 2))+
    geom_vline(aes(xintercept  = Start), data = filter(myMat_long,value == 1)) +
    coord_cartesian(ylim = c(0,1))+
    facet_wrap(~Name, ncol = 1, strip.position = 'left') +
    theme(axis.title.y       = element_blank(),
          axis.text.y        = element_blank(),
          axis.ticks.y       = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank())

pdf("BarcodePlot.pdf", width = 40, height = 20)
print(myChromosomePlot_plusStrand)
dev.off()

# add here colors for intergenic, gene, upstream or downstream!



# https://sib-swiss.github.io/enrichment-analysis-training/bonus_code/#code-for-a-lollipop-plot-with-ggplot2
# lollipop plot





