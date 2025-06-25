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

# read in the data
prX <- readRDS("SA6850_prXallWide_moreMetaInfo.rds")
mtX <- readRDS("methylation_data_EGGnogAnnotated.rds")

# here mtX centric
mtX_wPrX <- left_join(x = mtX, y = prX, by = c("locus_tag" = "LocusTag"))

# what catetories are there
table(mtX_wPrX$category)

# category of Interest
COI <- "gene"


# we want to filter out modified_base -> it is unclear what type of feature it is
mtX_wPrX <- mtX_wPrX[mtX_wPrX$feature != "modified_base",]
table(mtX_wPrX$feature)
table(mtX_wPrX$category)

# filter for Ancestor only
unique(mtX_wPrX$group)

# anyway we are only looking into ancestor (the rest is for students of Bio253)
mtX_wPrX <- mtX_wPrX[mtX_wPrX$group == "6850",]

tbl = mtX_wPrX %>% mutate(extracol = 1) %>% dplyr::select(1:5,111) %>% unique() %>% pivot_wider(names_from = "treatment", values_from = "extracol")
mtX_wPrX[mtX_wPrX$start %in% tbl[is.na(tbl$PASN) | is.na(tbl$TSB),]$start,] %>% group_by(group, treatment, feature, strand,category) %>% summarise(n_distinct(start)) %>% ggplot(aes(strand, `n_distinct(start)`, fill = treatment)) + geom_col(position = "dodge") + facet_grid(feature~category, scales = "free") + ggtitle("Sites unique to each treatment")

methylation_slim <- mtX_wPrX %>% select(treatment, start, strand, feature, locus_tag, category) |> distinct()
ggplot(methylation_slim, aes(x = start)) +
    geom_density(data = subset(methylation_slim, strand == "+"),
                 aes(y = after_stat(density), fill = strand),
                 alpha = 0.5, adjust = smoothin_param) +  # <-- less smoothing) +
    geom_density(data = subset(methylation_slim, strand == "-"),
                 aes(y = -after_stat(density), fill = strand),
                 alpha = 0.5,
                 adjust = smoothin_param) +  # <-- less smoothing) +
    scale_fill_manual(values = c("+" = "blue", "-" = "red")) +
    labs(
        title = "Density of Methylation Sites by Strand",
        x = "Chromosome Position",
        y = "Density",
        fill = "Strand"
    ) +
    theme_minimal() +
    geom_hline(yintercept = 0, linetype = "dashed") + facet_grid(category~paste0(treatment, feature))

methylation_slim %>% group_by(blocks, treatment, strand, feature,category) %>% summarise(n_distinct(start)) %>% pivot_wider(names_from = treatment, values_from = `n_distinct(start)`) %>% replace(is.na(.), 0) %>% mutate(diff = PASN - TSB) %>% mutate(ctt = case_when(diff > 0~ "hypomethylated in TSB",diff == 0 ~ "shared", diff < 0 ~"hypomethylated in PASN" )) %>% group_by(strand, feature, category, ctt) %>% summarise(number_of_blocks = n()) %>% ggplot(aes(strand, number_of_blocks, fill = ctt)) + geom_col(position = "dodge") + facet_grid(category ~ feature)
methylation_slim %>% group_by(blocks, treatment, strand, feature,category) %>% summarise(n_distinct(start)) %>% pivot_wider(names_from = treatment, values_from = `n_distinct(start)`) %>% replace(is.na(.), 0) %>% mutate(diff = PASN - TSB) %>% mutate(ctt = case_when(diff > 1~ "hypomethylated in TSB (>1 site)",.default ="shared", diff < -1 ~"hypomethylated in PASN (>1 site)" )) %>% group_by(strand, feature, category, ctt) %>% summarise(number_of_blocks = n()) %>% ggplot(aes(strand, number_of_blocks, fill = ctt)) + geom_col(position = "dodge") + facet_grid(category ~ feature) + ggtitle("Methylation density per 10kb")

