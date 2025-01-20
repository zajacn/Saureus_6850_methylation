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

# join in annotation from EggNogg and fastaCDS
# description lines from genomicCDS
desc <- read_delim("DescriptionLinesFromGenomicsFasta.txt", delim = "\t", col_names = FALSE)
str(desc)
# desc$proteinID <- gsub(".*protein_id=(.*?);.*", "\\1", desc$V1)
# ncbi description related
# xx$SAuniprotID <- sapply(strsplit(xx$SAorthologue, split = "_"), function(x)x[2])
# xx$locusTag <- gsub(x = xx$myMQprotGrps.firstDescription, pattern = ".*\\[locus_tag=(.*)\\] \\[protein=.*", replacement = "\\1")
# xx$proteinDesc <- gsub(x = xx$myMQprotGrps.firstDescription, pattern = ".*\\[protein=(.*)\\] \\[protein_id=.*", replacement = "\\1")
# xx$SAprotID <- gsub(x = xx$myMQprotGrps.firstDescription, pattern = ".*\\[protein_id=(.*)\\] \\[location=.*", replacement = "\\1")
# xx$SAlocation <- gsub(x = xx$myMQprotGrps.firstDescription, pattern = ".*\\[location=(.*)\\] \\[gbkey=.*", replacement = "\\1")

# DF
# initialize my df annoDF
desc$X1[1:101]
annoDF <- data.frame(matrix(ncol = 0, nrow = nrow(desc)))
annoDF$proteinID <- sapply(strsplit(gsub(x = sapply(strsplit(desc$X1, split = " "), function(x)x[1]), pattern = ">lcl\\|CP006706.1_cds_", replacement = ""),split = "_"), function(x)x[1])
annoDF$LocusTag <- gsub(x = desc$X1, pattern = ".*\\[locus_tag=(.*)\\] \\[protein=.*", replacement = "\\1")
annoDF$location <- gsub(x = desc$X1, pattern = ".*\\[location=(.*)\\] \\[gbkey=.*", replacement = "\\1")
#annoDF$proteinID <- gsub(x = desc$X1, pattern = ".*\\[protein_id=(.*)\\] \\[location=.*", replacement = "\\1")
annoDF$description <- gsub(x = desc$X1, , pattern = ".*\\[protein=(.*)\\] \\[protein_id=.*", replacement = "\\1")

# map in EGGnog results
# COG cat https://ecoliwiki.org/colipedia/index.php/Clusters_of_Orthologous_Groups_(COGs)
# http://eggnog-mapper.embl.de/MM_9y2o1_f1/
# mainly useful with out.emapper.annotations file

# load eggnog results
eggnog <- read_tsv("out.emapper.annotations", skip = 4)
str(eggnog)
# merge with annoDF
annoDF <- annoDF |> left_join(eggnog, by = c("proteinID" = "#query"))

# join main table
mtx <- mtx |> left_join(annoDF, by = c("locus_tag" = "LocusTag"))
