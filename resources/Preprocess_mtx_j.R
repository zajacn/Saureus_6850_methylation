#!/usr/local/bin/Rscript
#
# Jonas Grossmann <jg@fgcz.ethz.ch>
# 2025
#
#

# libs
library(dplyr)
library(readr)
library(tidyverse)


# load mtx
mtx <- readRDS("methylation_data.rds")

# do we see methylations in promotor regions? promotor regions are defined as 35bp downstream of the TSS
# https://www.researchgate.net/figure/Structure-of-bacterial-promoters-A-Schematic-of-RNA-polymerase-subunit-interactions_fig4_51866343

corePromotor <- c(-35, +20)
UpElement <- c(-65, -35)
PromotorWithUpelement <- c(-65, +20)

str(mtx)

# go strandwise
mtx |> filter(strand == "+") |> mutate(promotor = ifelse((start >= corePromotor[1] & start <= corePromotor[2]) | (end >= corePromotor[1] & end <= corePromotor[2]), 1, 0)) |> filter(promotor == 1)


# description lines from genomicCDS
desc <- read.table("DescriptionLinesFromGenomicsFasta.txt", sep = "\t", header = F, stringsAsFactors = F)
desc$V1[1]
# desc$proteinID <- gsub(".*protein_id=(.*?);.*", "\\1", desc$V1)
# ncbi description related
# xx$SAuniprotID <- sapply(strsplit(xx$SAorthologue, split = "_"), function(x)x[2])
# xx$locusTag <- gsub(x = xx$myMQprotGrps.firstDescription, pattern = ".*\\[locus_tag=(.*)\\] \\[protein=.*", replacement = "\\1")
# xx$proteinDesc <- gsub(x = xx$myMQprotGrps.firstDescription, pattern = ".*\\[protein=(.*)\\] \\[protein_id=.*", replacement = "\\1")
# xx$SAprotID <- gsub(x = xx$myMQprotGrps.firstDescription, pattern = ".*\\[protein_id=(.*)\\] \\[location=.*", replacement = "\\1")
# xx$SAlocation <- gsub(x = xx$myMQprotGrps.firstDescription, pattern = ".*\\[location=(.*)\\] \\[gbkey=.*", replacement = "\\1")



# map in EGGnog results
# COG cat https://ecoliwiki.org/colipedia/index.php/Clusters_of_Orthologous_Groups_(COGs)
# http://eggnog-mapper.embl.de/MM_9y2o1_f1/
# mainly useful with out.emapper.annotations file

# load eggnog results
eggnog <- read_tsv("out.emapper.annotations", skip = 4)

# merge
mtx <- mtx |> left_join(eggnog, by = c("gene_id" = "query_name"))
