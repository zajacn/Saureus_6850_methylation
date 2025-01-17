#!/usr/local/bin/Rscript
#
# Jonas Grossmann <jg@fgcz.ethz.ch>
# 2025
#
#

# load mtx

mtx <- readRDS("methylation_data.rds")

# do we see methylations in promotor regions? promotor regions are defined as 35bp downstream of the TSS
# https://www.researchgate.net/figure/Structure-of-bacterial-promoters-A-Schematic-of-RNA-polymerase-subunit-interactions_fig4_51866343

corePromotor <- c(-35, +20)
