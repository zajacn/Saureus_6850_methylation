#Load libraries
library(Biostrings)
library(DECIPHER)
library(stringdist)
library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(purrr)
library(rtracklayer)
library(ggseqlogo)
library(cluster)

# read in the data
prX <- readRDS("../resources/SA6850_prXallWide_moreMetaInfo.rds")
mtX <- readRDS("../resources/methylation_data_EGGnogAnnotated.rds")

# join the methylation and protein dataframes
mtX_wPrX <- left_join(x = mtX, y = prX, by = c("locus_tag" = "LocusTag"))

# we want to filter out modified_base -> it is unclear what type of feature it is
mtX_wPrX <- mtX_wPrX[mtX_wPrX$feature != "modified_base",]
mtX_wPrX <- mtX_wPrX[mtX_wPrX$feature != "m4C",]
#mtX_wPrX <- mtX_wPrX[mtX_wPrX$feature != "m6A",]

#Select motifs
motifs=mtX_wPrX$motif
# Compute pairwise string distances (Hamming works if all same length)
dist_mat <- stringdistmatrix(motifs, motifs, method = "hamming")

# Hierarchical clustering
clust <- hclust(as.dist(dist_mat), method = "average")

# Silhuette analysis
sil_widths <- sapply(100:120, function(k) {
    clusters <- cutree(clust, k = k)
    sil <- silhouette(clusters, as.dist(dist_mat))
    mean(sil[, 3])  # average silhouette width
})

plot(100:120, sil_widths, type = "b", pch = 19,
     xlab = "Number of clusters (k)",
     ylab = "Average silhouette width")

# Select the peak silhouette width
clusters <- cutree(clust, k = 113)

# Look at the clusters matrix
clusters

#Convert motifs to DNAStringSet
motif_set <- DNAStringSet(motifs)

# Derive consensus (IUPAC) for each cluster
cluster_cons <- tapply(1:length(motif_set), clusters, function(idx) {
    consensusString(motif_set[idx], ambiguityMap = IUPAC_CODE_MAP, threshold = 0.25)
})

# Have a look at the consesus sequences for each cluster
cluster_cons

# Select the top 10 motifs
top_scorers = data.frame(clusters) %>% group_by(clusters) %>% summarise(n()) %>% arrange(`n()`) %>% tail()

# Visualise logo plots for the top 10 motifs
for (i in top_scorers$clusters){
    print(ggseqlogo(motifs[clusters == i]))
}

# Lets merge the mtX prX table with the clusters to compare between groups and treatments
mtX_wPrX = cbind(mtX_wPrX, data.frame(clusters))

#Are the same motifs highly represented between treatments
mtX_wPrX %>% group_by(treatment, group, clusters) %>% summarise(n_distinct(start)) %>% ggplot(aes(treatment, `n_distinct(start)`, fill = as.character(clusters))) + geom_col(position = "dodge") + facet_grid(~group)

#Is there a difference in motif usage between categories
mtX_wPrX %>% group_by(treatment, group, category, clusters) %>% summarise(n_distinct(start)) %>% ggplot(aes(treatment, `n_distinct(start)`, fill = as.character(clusters))) + geom_col(position = "dodge") + facet_grid(category~group)

