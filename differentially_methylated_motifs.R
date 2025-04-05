#This R script is a collection of code for looking into motifs methylated in TSB or in PASN but not vice versa

mtx = readRDS("/home/zajac/Saureus_6850_methylation/resources/methylation_data.rds")
missingTSB = mtx[mtx$feature == "m4C",] %>% 
  group_by(group, treatment, motif) %>% 
  dplyr::summarise(motifs = n_distinct(start)) %>% 
  filter(motifs > 1) %>% 
  pivot_wider(names_from = treatment, values_from = motifs) %>% 
  filter(is.na(TSB))
missingTSB = missingTSB$motif
missingPASN = mtx[mtx$feature == "m4C",] %>% 
  group_by(group, treatment, motif) %>% 
  dplyr::summarise(motifs = n_distinct(start)) %>% 
  filter(motifs > 1) %>% 
  pivot_wider(names_from = treatment, values_from = motifs) %>% 
  filter(is.na(PASN))
missingPASN = missingPASN$motif

#Genes with motifs methylated in PASN but not in TSB
unique(mtx[mtx$motif %in% missingTSB,]$ID)

#Genes with motifs methylated in TSB but not in PASN
unique(mtx[mtx$motif %in% missingTSB,]$ID)